/* sim.js -- browser/Node port of the V3 solvers used by the interactives.
   Mirrors code/params.py, code/twodomain.py and code/v3model.py (FisherKPP).
   Units: cm, s, mol/cm3 unless a name says otherwise. */
var Sim = (function () {
  "use strict";
  var UM = 1e-9, DAY = 86400, HOUR = 3600, LN2 = Math.LN2;

  var P = {
    D0_V14: 4.5e-6, A_V14: 0.73, A_FIBRE: 0.6, XI: 14, EPS_GEL: 0.99,
    EPS_TIS: 0.8, D_TIS: 1.2e-6, L_GEL: 0.05, L_TIS: 0.15, K_P: 1,
    T_HALF_V14: 6.1 * HOUR, C_SOL_UM: 738, TARGET_UM: 25, DELTA: 0.05,
    // closure (Model 6)
    D_N: 1e-9, R0: 0.8 / DAY, D0: 0.24 / DAY, DD: 0.3 / DAY, E_MAX: 2, EC50: 0.5,
    BETA: 0.5, K_IR: 0.3, K_ID: 0.3, S_I: 2, KE0: 1 / (18 * HOUR), W: 75e-4,
    KD_L_UM: 0.065, L_UM: 1
  };
  P.K_PROT = LN2 / P.T_HALF_V14;
  P.K_DEG = LN2 / P.T_HALF_V14;

  function amsdenH(a, xi, af) {
    var q = (a + af) / (xi + 2 * af);
    return Math.exp(-Math.PI * q * q);
  }
  function gelD(o) {
    o = o || {};
    var D0 = o.D0 != null ? o.D0 : P.D0_V14, eps = o.eps != null ? o.eps : P.EPS_GEL;
    return D0 * amsdenH(P.A_V14, P.XI, P.A_FIBRE) / Math.pow(eps, -0.5);
  }

  /* ---------- tridiagonal helpers ---------- */
  function Tri(sub, dia, sup, dt) {
    var n = dia.length, cp = new Float64Array(n), inv = new Float64Array(n);
    this.n = n; this.lo = new Float64Array(n); this.cp = cp; this.inv = inv;
    for (var i = 0; i < n; i++) this.lo[i] = -dt * sub[i];
    var up = new Float64Array(n);
    for (i = 0; i < n; i++) up[i] = -dt * sup[i];
    var d = 1 - dt * dia[0];
    inv[0] = 1 / d; cp[0] = up[0] * inv[0];
    for (i = 1; i < n; i++) {
      d = 1 - dt * dia[i] - this.lo[i] * cp[i - 1];
      inv[i] = 1 / d; cp[i] = up[i] * inv[i];
    }
  }
  Tri.prototype.solve = function (b, x) {
    var n = this.n, lo = this.lo, cp = this.cp, inv = this.inv, i;
    x[0] = b[0] * inv[0];
    for (i = 1; i < n; i++) x[i] = (b[i] - lo[i] * x[i - 1]) * inv[i];
    for (i = n - 2; i >= 0; i--) x[i] -= cp[i] * x[i + 1];
  };

  /* ---------- Models 1-2: monolithic gel + tissue ---------- */
  var SEG = [[0.5, 120], [5, 108], [30, 2800], [300, 600]];

  function twoDomain(o) {
    o = o || {};
    var Lgel = o.Lgel != null ? o.Lgel : P.L_GEL, Ltis = o.Ltis != null ? o.Ltis : P.L_TIS;
    var RV = o.RV != null ? o.RV : 1, kprot = o.kprot != null ? o.kprot : P.K_PROT;
    var Dg = o.Dg != null ? o.Dg : gelD(), Dt = o.Dt != null ? o.Dt : P.D_TIS;
    var CVtot = (o.CVtotUM != null ? o.CVtotUM : P.C_SOL_UM) * UM;
    var nxg = o.nxg || 40, nxt = o.nxt || 120, tEnd = (o.tEndDays || 3) * DAY;
    var delta = o.delta != null ? o.delta : P.DELTA, target = (o.targetUM || P.TARGET_UM) * UM;
    var Kp = P.K_P, eg = P.EPS_GEL, et = P.EPS_TIS, kdeg = P.K_DEG;
    var hg = Lgel / nxg, ht = Ltis / nxt, ng = nxg, nt = nxt - 1, N = ng + nt;
    var a = eg * Dg / hg, b = et * Dt / ht, den = b + a * Kp, ca = a / den, cb = b / den;
    var sub = new Float64Array(N), dia = new Float64Array(N), sup = new Float64Array(N), i, j, p;
    var ct = Dt / (ht * ht), cg = Dg / (hg * hg);
    for (i = 1; i <= nt; i++) {                     // tissue node t_i at pos nt-i
      p = nt - i;
      dia[p] = -2 * ct - kprot;
      if (i < nt) sub[p] = ct;
      if (i > 1) sup[p] = ct; else { sup[p] = ct * ca; dia[p] += ct * cb; }
    }
    for (j = 1; j <= ng; j++) {                     // gel node g_j at pos nt+j-1
      p = nt + j - 1;
      var d0 = -2 * cg - kdeg, lo = cg, up = cg;
      if (j === 1) { d0 += cg * Kp * ca; lo = cg * Kp * cb; }
      if (j === ng) { up = 0; if (j > 1) lo = 2 * cg; }
      dia[p] = d0 / RV; sub[p] = lo / RV; sup[p] = up / RV;
    }
    var C = new Float64Array(N), X = new Float64Array(N);
    var cfree = CVtot / RV;
    for (j = 1; j <= ng; j++) C[nt + j - 1] = cfree;
    var kmax = Math.floor(delta / ht + 1e-9), ts = [0], avg = [0], t = 0, done = false;
    var peak = 0, peakT = 0;
    function depthAvg() {
      var c0 = ca * C[nt] + cb * C[nt - 1], s = 0.5 * c0;
      for (var k = 1; k < kmax; k++) s += C[nt - k];
      s += 0.5 * (kmax <= nt ? C[nt - kmax] : 0);
      return s * ht / delta;
    }
    for (var sI = 0; sI < SEG.length && !done; sI++) {
      var dt = SEG[sI][0], tri = new Tri(sub, dia, sup, dt);
      for (var s = 0; s < SEG[sI][1]; s++) {
        tri.solve(C, X); var tmp = C; C = X; X = tmp; t += dt;
        var v = depthAvg(); ts.push(t); avg.push(v);
        if (v > peak) { peak = v; peakT = t; }
        if (t >= tEnd) { done = true; break; }
      }
    }
    var dur = 0;
    for (i = 1; i < ts.length; i++) if (avg[i] >= target || avg[i - 1] >= target) dur += ts[i] - ts[i - 1];
    return {
      t: ts, avg: avg, hoursAbove: dur / HOUR, peakUM: peak / UM, peakT: peakT,
      totalUM: CVtot / UM, freeUM: cfree / UM, solFrac: CVtot / UM / P.C_SOL_UM,
      overCeiling: CVtot / UM > P.C_SOL_UM + 1e-9
    };
  }

  /* ---------- closure rates from the register closures ---------- */
  function rpOf(Ce, Ie) {
    return P.R0 * (1 + P.E_MAX * Ce / (P.EC50 + Ce)) * (1 - P.BETA * Ie / (P.K_IR + Ie));
  }
  function dOf(Ie) {
    var q = Math.pow(Ie, P.S_I);
    return P.D0 + P.DD * q / (Math.pow(P.K_ID, P.S_I) + q);
  }
  var a_lps = P.L_UM / P.KD_L_UM, theta0 = a_lps / (1 + a_lps);
  var RATES = {
    rpChr: rpOf(0, theta0), dChr: dOf(theta0),
    rpTrt: rpOf(P.EC50, theta0 / 2), dTrt: dOf(theta0 / 2)
  };
  RATES.thr = 0.5 * (1 - RATES.dTrt / RATES.rpTrt);   // pinned front level

  function frontAt(n, thr, r, dr, nr) {
    if (n[0] >= thr) return 0;
    for (var i = 1; i <= nr; i++) if (n[i] >= thr) {
      var n0 = n[i - 1], n1 = n[i];
      if (n1 === n0) return r[i];
      return r[i - 1] + (thr - n0) / (n1 - n0) * dr;
    }
    return r[nr];
  }

  /* Radial Fisher-KPP with a death term; sFn(tSeconds) in [0,1] drives the
     drug effect.  opts.frames > 0 stores snapshots for the animation. */
  function fisher(o) {
    var Rw = o.Rw, nr = o.nr || 300, Rinf = Rw * 2.2, dr = Rinf / nr, K = 1;
    var tEnd = o.tEndDays * DAY, dt = 0.2 * dr * dr / P.D_N, thr = RATES.thr;
    var n = new Float64Array(nr + 1), m = new Float64Array(nr + 1), r = new Float64Array(nr + 1), i;
    for (i = 0; i <= nr; i++) { r[i] = i * dr; n[i] = 0.5 * (1 + Math.tanh((r[i] - Rw) / P.W)); }
    var steps = Math.floor(tEnd / dt), D = P.D_N, rec = 40,
        frames = [], ts = [], fr = [], closed = false, tc = Infinity, t = 0, frameEvery = o.frameEvery ? o.frameEvery * DAY : 0, nextF = 0;
    var sGrid = o.sGrid, sdt = o.sdt || 60;
    for (var st = 1; st <= steps; st++) {
      var s = sGrid ? sGrid[Math.min(Math.floor(t / sdt), sGrid.length - 1)] : o.sFn(t);
      var rp = RATES.rpChr + (RATES.rpTrt - RATES.rpChr) * s, dd = RATES.dChr + (RATES.dTrt - RATES.dChr) * s;
      var c0 = D / (dr * dr);
      m[0] = n[0] + dt * (c0 * 4 * (n[1] - n[0]) + rp * n[0] * (1 - n[0] / K) - dd * n[0]);
      for (i = 1; i < nr; i++) {
        var lap = c0 * (n[i + 1] - 2 * n[i] + n[i - 1]) + D * (n[i + 1] - n[i - 1]) / (2 * dr * r[i]);
        m[i] = n[i] + dt * (lap + rp * n[i] * (1 - n[i] / K) - dd * n[i]);
      }
      for (i = 0; i < nr; i++) n[i] = m[i] > 0 ? m[i] : 0;
      n[nr] = K; t += dt;
      if (st % rec === 0 || st === steps) {
        var rf = frontAt(n, thr, r, dr, nr);
        ts.push(t); fr.push(rf);
        if (rf <= 1e-9) { closed = true; tc = t; if (!o.keepGoing) break; }
      }
      if (frameEvery && t >= nextF) { frames.push(Float32Array.from(n)); nextF += frameEvery; }
    }
    return { closed: closed, tc: tc, t: ts, rf: fr, r: r, frames: frames, Rinf: Rinf, finalN: n };
  }

  /* effect state s(t) on a 60 s grid for a patch changed every `intervalH` h
     that keeps tissue above target for `covH` h */
  function scheduleS(intervalH, covH, tEndDays) {
    var dt = 60, nn = Math.ceil(tEndDays * DAY / dt) + 2, s = new Float64Array(nn), e = Math.exp(-dt * P.KE0);
    for (var i = 1; i < nn; i++) {
      var u = (((i - 1) * dt / HOUR) % intervalH) < covH ? 1 : 0;
      s[i] = u + (s[i - 1] - u) * e;
    }
    return s;
  }
  function closureTimeInterval(Rw, intervalH, covH, tEndDays, nr) {
    var T = tEndDays || 120;
    var res = fisher({ Rw: Rw, nr: nr || 300, tEndDays: T, sGrid: scheduleS(intervalH, covH, T) });
    return res.closed ? res.tc / DAY : Infinity;
  }
  function maxInterval(Rw, covH, hi, nr) {
    hi = hi || 24;
    if (!(covH > 0)) return 0;
    if (isFinite(closureTimeInterval(Rw, hi, covH, 120, nr))) return hi;
    var lo = covH;
    if (!isFinite(closureTimeInterval(Rw, lo * 1.0001, covH, 120, nr))) return 0;
    while (hi - lo > 0.25) {
      var mid = 0.5 * (lo + hi);
      if (isFinite(closureTimeInterval(Rw, mid, covH, 120, nr))) lo = mid; else hi = mid;
    }
    return lo;
  }
  /* Same bisection as maxInterval, but one closure solve per task (~20 ms
     each) so the browser can keep scrolling between steps. Returns a cancel
     function; done(hours) is called once with the result. */
  function maxIntervalAsync(Rw, covH, done, hi, nr) {
    hi = hi || 24; nr = nr || 300;
    var cancelled = false, lo = covH, stage = 0;
    function ok(h) { return isFinite(closureTimeInterval(Rw, h, covH, 120, nr)); }
    function later(f) { setTimeout(function () { if (!cancelled) f(); }, 0); }
    function step() {
      if (!(covH > 0)) return done(0);
      if (stage === 0) { stage = 1; if (ok(hi)) return done(hi); return later(step); }
      if (stage === 1) { stage = 2; if (!ok(lo * 1.0001)) return done(0); return later(step); }
      if (hi - lo <= 0.25) return done(lo);
      var mid = 0.5 * (lo + hi);
      if (ok(mid)) lo = mid; else hi = mid;
      later(step);
    }
    later(step);
    return function () { cancelled = true; };
  }
  function pulseS(Tdays) {
    var T = Tdays * DAY, tau = 1 / P.KE0;
    return function (t) { return t <= T ? 1 : Math.exp(-(t - T) / tau); };
  }

  /* ---------- Models 3-4 closed forms ---------- */
  function suppression(psi, a) { return (1 + a) / (1 + a + psi); }
  function thetaLps(psi, a) { return a / (1 + a + psi); }
  function p65(psi, a, gamma, h) {
    var x = Math.pow(a * gamma, h);
    return x / (x + Math.pow(1 + a + psi, h));
  }

  return {
    P: P, UM: UM, DAY: DAY, HOUR: HOUR, RATES: RATES, gelD: gelD, twoDomain: twoDomain,
    fisher: fisher, frontAt: frontAt, scheduleS: scheduleS, closureTimeInterval: closureTimeInterval,
    maxInterval: maxInterval, maxIntervalAsync: maxIntervalAsync, pulseS: pulseS, suppression: suppression,
    thetaLps: thetaLps, p65: p65
  };
})();
if (typeof module !== "undefined") module.exports = Sim;
