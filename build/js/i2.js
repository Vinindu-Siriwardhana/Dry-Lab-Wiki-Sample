/* Interactive 2 -- closure animation: radial Fisher-KPP, n(r,t).
   Three modes: repeated patch changes (each patch covers COV h, the effect
   compartment fades in between), one long pulse, or untreated.

   The wound bed is drawn as tissue rather than as a colour ramp. Every
   fibroblast mark sits at a fixed site with its own threshold u in (0,1) and is
   shown only where the local density n(r)/K exceeds u, so the fraction of marks
   present at any radius IS the model density -- the picture stays a plot. */
(function () {
  "use strict";
  function el(id) { return document.getElementById(id); }
  var disc = el("i2-disc");
  if (!disc) return;
  var S = Sim, C = Plot.COL, prof = el("i2-prof"), tl = el("i2-time");
  var DAYMAX = 90, FE = 0.25, SDT = 60, NR = 300, TAU = 6.2832;
  var K = {
    bed: [194, 72, 79], dermis: [239, 205, 189],          // open granulation bed -> repaired tissue
    cell: "#FBEBD6", cellEdge: "rgba(150,98,62,.55)", nucleus: "rgba(150,98,62,.75)",
    cellLine: "#9C6B45", cellFill: "rgba(185,138,99,.16)",
    front: "#16201B", edge: "#8C6A57", drug: "33,102,172", thr: C.thresh
  };
  var st = { mode: "patch", T: 14, iv: 24, Rw: 0.2, res: null, rf: [], sArr: null, band: null,
             endDay: DAYMAX, xmax: DAYMAX, day: 0, playing: false, last: 0, visible: false,
             parked: false, autoplayed: false, speed: "normal", flashAt: -1e9, lastPatch: -1 };
  var COV = null;   // hours one patch keeps tissue above target (500 um alginate, defaults of Interactive 3)
  var reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  var nBuf = new Float32Array(NR + 1);

  function mulberry32(a) {
    return function () {
      a |= 0; a = a + 0x6D2B79F5 | 0;
      var t = Math.imul(a ^ a >>> 15, 1 | a);
      t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }
  function clamp(x, a, b) { return x < a ? a : x > b ? b : x; }
  function smooth(x) { x = clamp(x, 0, 1); return x * x * (3 - 2 * x); }

  /* Cell sites in unit-disc coordinates. The outer domain is always
     R_inf = 2.2 R_w, so the normalised geometry -- and therefore this set -- is
     the same for every wound radius. Cells inside the original wound point
     radially, the direction they migrate; cells in intact tissue lie at random. */
  var NO = 16, SZ = [0.88, 1, 1.12];   // sprite orientations over [0, pi) and size classes
  var SITES = (function () {
    var rnd = mulberry32(20261015), out = [], g = 0.046, x, y;
    for (y = -1; y <= 1; y += g) for (x = -1; x <= 1; x += g) {
      var px = x + (rnd() - 0.5) * g * 0.9, py = y + (rnd() - 0.5) * g * 0.9, rr = Math.sqrt(px * px + py * py);
      if (rr > 0.975) continue;
      var radial = rr < 0.5;
      var k = 0.85 + 0.3 * rnd(), a = radial ? Math.atan2(py, px) + (rnd() - 0.5) * 0.7 : rnd() * Math.PI;
      out.push({ x: px, y: py, r: rr, u: 0.03 + 0.94 * rnd(),
                 z: k < 0.95 ? 0 : k < 1.05 ? 1 : 2,
                 o: Math.round(((a % Math.PI) + Math.PI) % Math.PI / Math.PI * NO) % NO });
    }
    return out;
  })();

  function coverage() {
    if (COV === null) COV = S.twoDomain({ CVtotUM: S.P.C_LOAD_V14_UM, Lgel: S.P.L_GEL, RV: 1, kprot: S.P.K_PROT }).hoursAbove;
    return COV;
  }
  function niceMax(d) {
    var steps = [10, 15, 20, 25, 30, 40, 50, 60, 75, 90];
    for (var i = 0; i < steps.length; i++) if (steps[i] >= d) return steps[i];
    return DAYMAX;
  }

  /* ---------- simulation ---------- */
  function simulate() {
    var o = { Rw: st.Rw, nr: NR, tEndDays: DAYMAX, frameEvery: FE }, nS = Math.ceil(DAYMAX * S.DAY / SDT) + 2, i;
    if (st.mode === "patch") { st.sArr = S.scheduleS(st.iv, coverage(), DAYMAX); o.sGrid = st.sArr; o.sdt = SDT; }
    else {
      var f = st.mode === "pulse" ? S.pulseS(st.T) : function () { return 0; };
      o.sFn = f; st.sArr = new Float64Array(nS);
      for (i = 0; i < nS; i++) st.sArr[i] = f(i * SDT);
    }
    var res = S.fisher(o), dr = res.Rinf / NR, thr = S.RATES.thr;
    st.res = res;
    st.rf = res.frames.map(function (n) { return S.frontAt(n, thr, res.r, dr, NR); });
    st.endDay = res.closed ? res.tc / S.DAY : (res.frames.length - 1) * FE;
    st.xmax = res.closed ? niceMax(st.endDay * 1.12) : DAYMAX;
    st.band = buildBand();
    el("i2-d").max = Math.floor(st.endDay * 4) / 4;
    st.day = 0; st.lastPatch = -1;
  }

  /* drug effect s(t) binned to ~2 px columns: min/max envelope of the
     saw-tooth plus its trailing 1-day mean, which is the quantity the 33 %
     closure condition is about */
  function buildBand() {
    var M = 360, s = st.sArr, n = s.length, cum = new Float64Array(n + 1), i, j, W = Math.round(S.DAY / SDT);
    for (i = 0; i < n; i++) cum[i + 1] = cum[i] + s[i];
    var xs = [], lo = [], hi = [], mean = [];
    for (i = 0; i < M; i++) {
      var d0 = i / M * st.xmax, d1 = (i + 1) / M * st.xmax;
      var a = Math.min(n - 1, Math.floor(d0 * S.DAY / SDT)), b = Math.min(n, Math.max(a + 1, Math.ceil(d1 * S.DAY / SDT)));
      var mn = 1, mx = 0;
      for (j = a; j < b; j++) { if (s[j] < mn) mn = s[j]; if (s[j] > mx) mx = s[j]; }
      var e = b, w0 = Math.max(0, e - W);
      xs.push(0.5 * (d0 + d1)); lo.push(mn); hi.push(mx); mean.push((cum[e] - cum[w0]) / Math.max(1, e - w0));
    }
    return { x: xs, lo: lo, hi: hi, mean: mean };
  }
  function sAt(day) { return st.sArr[Math.min(st.sArr.length - 1, Math.max(0, Math.floor(day * S.DAY / SDT)))]; }
  function nssAt(s) {
    var R = S.RATES, rp = R.rpChr + (R.rpTrt - R.rpChr) * s, d = R.dChr + (R.dTrt - R.dChr) * s;
    return Math.max(0, 1 - d / rp);
  }

  /* density and front interpolated between the stored 0.25-day frames, so
     the animation moves continuously instead of stepping */
  function current() {
    var fr = st.res.frames, last = fr.length - 1, x = clamp(st.day / FE, 0, last),
        k0 = Math.floor(x), k1 = Math.min(last, k0 + 1), w = x - k0, a = fr[k0], b = fr[k1], i;
    for (i = 0; i <= NR; i++) nBuf[i] = a[i] + (b[i] - a[i]) * w;
    return { n: nBuf, rf: st.rf[k0] + (st.rf[k1] - st.rf[k0]) * w, k: k0 };
  }
  function patchesBy(day) { return Math.floor(day * 24 / st.iv) + 1; }
  function speedDays() {
    var normal = Math.max(1.5, st.endDay / 12);
    return st.speed === "slow" ? 0.25 : st.speed === "fast" ? normal * 3 : normal;
  }

  /* ---------- drawing ---------- */
  function draw(now) {
    if (!st.res) return;
    var cur = current();
    el("i2-d-v").textContent = st.day.toFixed(1) + " d";
    drawDisc(cur, now || 0);
    drawProfile(cur);
    drawTimeline(cur);
    status(cur);
  }

  /* One fibroblast pre-rendered per (size, orientation): blitting ~1500 small
     sprites costs a fraction of filling and stroking ~1500 ellipse paths, and
     the cell field is the only part of the disc redrawn every frame. */
  var spr = null, sprKey = "";
  function sprites(a0, dpr) {
    var key = a0.toFixed(3) + "|" + dpr;
    if (spr && sprKey === key) return spr;
    var out = [], zi, o;
    for (zi = 0; zi < SZ.length; zi++) {
      var row = [];
      for (o = 0; o < NO; o++) {
        var ra = a0 * SZ[zi], rb = ra / 2.6, ang = o / NO * Math.PI, half = ra + 2,
            px = Math.ceil(2 * half * dpr), cv = document.createElement("canvas"), c;
        cv.width = cv.height = px; c = cv.getContext("2d");
        c.setTransform(dpr, 0, 0, dpr, px / 2, px / 2);
        c.beginPath(); c.ellipse(0, 0, ra, rb, ang, 0, TAU); c.fillStyle = K.cell; c.fill();
        c.lineWidth = 0.6; c.strokeStyle = K.cellEdge; c.stroke();
        c.beginPath(); c.ellipse(0, 0, ra / 3.4, ra / 3.4 / 1.7, ang, 0, TAU); c.fillStyle = K.nucleus; c.fill();
        row.push({ cv: cv, half: half });
      }
      out.push(row);
    }
    spr = out; sprKey = key; return spr;
  }

  /* The cell field is rendered into its own layer at native device size and
     refreshed at most every 50 ms while playing (~20 Hz). Cells fade in over
     several frames, so the cap is invisible, while the front, rings and haze
     still move every frame. Scrubbing and static redraws always refresh it. */
  var layer = null, layerAt = -1e9;
  function cellLayer(n, cx, cy, R, dpr, bw, bh, now) {
    var fresh = !layer || layer.width !== bw || layer.height !== bh;
    if (fresh) { layer = document.createElement("canvas"); layer.width = bw; layer.height = bh; }
    if (!fresh && now && now - layerAt < 50) return layer;
    var c = layer.getContext("2d"), SP = sprites(0.052 * R, dpr), i;
    c.setTransform(1, 0, 0, 1, 0, 0); c.clearRect(0, 0, bw, bh);
    for (i = 0; i < SITES.length; i++) {
      var p = SITES[i], vis = smooth((n[Math.min(NR, Math.round(p.r * NR))] - p.u) / 0.07 + 0.5);
      if (vis <= 0.08) continue;
      var sp = SP[p.z][p.o], w = sp.cv.width;
      c.globalAlpha = vis;     // whole device pixels at native size: the rasteriser's fast copy path
      c.drawImage(sp.cv, Math.round((cx + p.x * R) * dpr - w / 2), Math.round((cy + p.y * R) * dpr - w / 2));
    }
    c.globalAlpha = 1;
    layerAt = now || -1e9;
    return layer;
  }

  var tex = null, texKey = "";
  function texture(px) {
    if (tex && texKey === px) return tex;
    var cv = document.createElement("canvas"), rnd = mulberry32(7), c, i, n = Math.round(px * px / 46);
    cv.width = cv.height = px; c = cv.getContext("2d");
    for (i = 0; i < n; i++) {
      c.fillStyle = rnd() < 0.55 ? "rgba(90,18,28,.16)" : "rgba(255,238,228,.14)";
      c.beginPath(); c.arc(rnd() * px, rnd() * px, 0.5 + rnd() * 1.5 * px / 600, 0, TAU); c.fill();
    }
    tex = cv; texKey = px; return tex;
  }
  function rgb(a, b, t) {
    return "rgb(" + Math.round(a[0] + (b[0] - a[0]) * t) + "," + Math.round(a[1] + (b[1] - a[1]) * t) + "," +
      Math.round(a[2] + (b[2] - a[2]) * t) + ")";
  }
  function pill(c, x, y, txt, bg, fg) {
    c.font = "600 11px Inter, system-ui, sans-serif";
    var w = c.measureText(txt).width + 14, h = 20, r = 10;
    c.beginPath(); c.moveTo(x + r, y); c.lineTo(x + w - r, y); c.arc(x + w - r, y + r, r, -Math.PI / 2, Math.PI / 2);
    c.lineTo(x + r, y + h); c.arc(x + r, y + r, r, Math.PI / 2, 1.5 * Math.PI); c.closePath();
    c.fillStyle = bg; c.fill();
    c.fillStyle = fg; c.textAlign = "left"; c.textBaseline = "middle"; c.fillText(txt, x + 7, y + h / 2 + 0.5);
    return w;
  }
  function ring(c, cx, cy, rad, col, w, dash) {
    c.beginPath(); c.strokeStyle = col; c.lineWidth = w; c.setLineDash(dash || []);
    c.arc(cx, cy, Math.max(0.1, rad), 0, TAU); c.stroke(); c.setLineDash([]);
  }

  function drawDisc(cur, now) {
    var s = Plot.setup(disc, 1), c = s.c, cx = s.w / 2, cy = s.h / 2, R = s.w / 2 - 21, i,
        n = cur.n, Rinf = st.res.Rinf, dpr = window.devicePixelRatio || 1, sNow = sAt(st.day);

    // tissue: open bed blends into repaired tissue with the local density
    var g = c.createRadialGradient(cx, cy, 0, cx, cy, R), STOPS = 64;
    for (i = 0; i <= STOPS; i++) g.addColorStop(i / STOPS, rgb(K.bed, K.dermis, smooth((n[Math.round(i / STOPS * NR)] - 0.04) / 0.5)));
    c.save(); c.beginPath(); c.arc(cx, cy, R, 0, TAU); c.clip();
    c.fillStyle = g; c.fillRect(cx - R, cy - R, 2 * R, 2 * R);
    c.drawImage(texture(Math.round(2 * R * dpr)), cx - R, cy - R, 2 * R, 2 * R);

    // fibroblasts: a site is shown where the local density exceeds its threshold
    var lay = cellLayer(n, cx, cy, R, dpr, disc.width, disc.height, now);
    c.save(); c.setTransform(1, 0, 0, 1, 0, 0); c.drawImage(lay, 0, 0); c.restore();

    // drug effect: a blue haze whose strength is the effect compartment s(t)
    if (sNow > 0.004) { c.fillStyle = "rgba(" + K.drug + "," + (0.2 * sNow).toFixed(3) + ")"; c.fillRect(cx - R, cy - R, 2 * R, 2 * R); }
    c.restore();

    // patch halo; a flash marks each new patch when playback is slow enough to see one
    if (sNow > 0.004) ring(c, cx, cy, R + 5, "rgba(" + K.drug + "," + (0.25 + 0.6 * sNow).toFixed(3) + ")", 4);
    var fl = reduce ? 0 : clamp(1 - (now - st.flashAt) / 520, 0, 1);
    if (fl > 0) ring(c, cx, cy, R + 5 + 9 * (1 - fl), "rgba(" + K.drug + "," + (0.7 * fl).toFixed(3) + ")", 3 + 5 * fl);
    else ring(c, cx, cy, R + 1, "#D5DBE0", 1.5);

    // original wound edge and the healing front
    var re = st.Rw / Rinf * R;
    ring(c, cx, cy, re, "rgba(255,255,255,.85)", 3, [5, 4]);
    ring(c, cx, cy, re, K.edge, 1.3, [5, 4]);
    if (cur.rf > 1e-6) { ring(c, cx, cy, cur.rf / Rinf * R, "rgba(255,255,255,.95)", 5.5); ring(c, cx, cy, cur.rf / Rinf * R, K.front, 2.4); }

    // scale bar: 1 mm
    var bar = 0.1 / Rinf * R, bx = 10, by = s.h - 10;
    c.strokeStyle = K.front; c.lineWidth = 2.5; c.beginPath(); c.moveTo(bx, by); c.lineTo(bx + bar, by); c.stroke();
    c.font = "600 11px Inter, system-ui, sans-serif"; c.fillStyle = K.front; c.textAlign = "center"; c.textBaseline = "bottom";
    c.fillText("1 mm", bx + bar / 2, by - 4);

    // labels
    pill(c, 6, 6, "Day " + st.day.toFixed(1), "rgba(255,255,255,.95)", "#16201B");
    var tag = st.mode === "patch" ? ["Patch #" + patchesBy(st.day), "rgba(" + K.drug + ",.95)"]
      : st.mode === "pulse" ? (st.day <= st.T ? ["Treating", "rgba(" + K.drug + ",.95)"] : ["Treatment over", "rgba(107,116,125,.95)"])
      : ["No treatment", "rgba(107,116,125,.95)"];
    c.font = "600 11px Inter, system-ui, sans-serif";
    pill(c, s.w - c.measureText(tag[0]).width - 20, 6, tag[0], tag[1], "#fff");
    if (st.res.closed && st.day >= st.endDay - 1e-6) {
      c.font = "600 11px Inter, system-ui, sans-serif";
      pill(c, s.w - c.measureText("Closed").width - 20, s.h - 26, "Closed", "#1B6B45", "#fff");
    }
  }

  function drawProfile(cur) {
    var s = Plot.setup(prof, 1), n = cur.n, r = st.res.r, mmax = st.res.Rinf * 10, xs = [], ys = [], i,
        step = Math.max(1, Math.floor(r.length / 160)), tk = [], thr = S.RATES.thr, nss = nssAt(sAt(st.day));
    for (var v = 0; v <= mmax; v += mmax > 8 ? 2 : 1) tk.push([v, String(v)]);
    var ax = Plot.Axes(s, { l: 44, t: 14, r: 12, b: 48 }, [0, mmax], [0, 1.05], {
      xlabel: "distance from wound centre (mm)", ylabel: "cell density n / K", xticks: tk,
      yticks: [[0, "0"], [0.5, "0.5"], [1, "1"]]
    });
    ax.frame();
    for (i = 0; i < r.length; i += step) { xs.push(r[i] * 10); ys.push(n[i]); }
    xs.push(r[r.length - 1] * 10); ys.push(n[n.length - 1]);
    var c = ax.c;
    c.beginPath(); c.moveTo(ax.X(xs[0]), ax.Y(0));
    for (i = 0; i < xs.length; i++) c.lineTo(ax.X(xs[i]), ax.Y(ys[i]));
    c.lineTo(ax.X(xs[xs.length - 1]), ax.Y(0)); c.closePath(); c.fillStyle = K.cellFill; c.fill();
    ax.line(xs, ys, K.cellLine, 2.6);
    ax.vline(st.Rw * 10, K.edge, [4, 4], 1.2);

    // what the current drug level can sustain, against what counts as closed
    ax.hline(thr, K.thr, [5, 4], 1.4);
    ax.hline(Math.max(0.006, nss), "rgb(" + K.drug + ")", [2, 3], 1.8);
    var f = "600 10.5px Inter, system-ui, sans-serif", kx = ax.x0 + 8, ky = ax.y1 + 8;
    function key(y, col, dash, w, txt) {
      c.beginPath(); c.strokeStyle = col; c.lineWidth = w; c.setLineDash(dash);
      c.moveTo(kx, y); c.lineTo(kx + 20, y); c.stroke(); c.setLineDash([]);
      ax.text(txt, kx + 26, y, col, "left", f);
    }
    c.font = f; var k1 = nss > 0.006 ? "max at current drug level" : "drug too low: cells die back";
    c.fillStyle = "rgba(255,255,255,.88)";
    c.fillRect(kx - 4, ky - 9, 34 + Math.max(c.measureText(k1).width, c.measureText("closure threshold").width), 34);
    key(ky, "rgb(" + K.drug + ")", [2, 3], 1.8, k1);
    key(ky + 16, K.thr, [5, 4], 1.4, "closure threshold");
    if (cur.rf > 1e-6) ax.dot(cur.rf * 10, thr, K.front, 4.5);
  }

  function drawTimeline(cur) {
    var s = Plot.setup(tl, 0.42), c = s.c, xmax = st.xmax, b = st.band, i;
    var T = 10, Bm = 38, gap = 16, laneH = Math.round((s.h - T - Bm - gap) * 0.36);
    var xt = [], step = xmax <= 20 ? 5 : xmax <= 50 ? 10 : 15;
    for (var d = 0; d <= xmax + 1e-9; d += step) xt.push([d, String(d)]);
    var play = Math.min(st.day, st.endDay), lbl = "600 10.5px Inter, system-ui, sans-serif";

    /* lane 1: drug effect */
    var A = Plot.Axes(s, { l: 50, t: T, r: 14, b: s.h - T - laneH }, [0, xmax], [0, 1.05], {
      ylabel: "drug", xticks: xt.map(function (t) { return [t[0], ""]; }), yticks: [[0, "0"], [1, "max"]]
    });
    A.frame();
    function band(alpha, lineAlpha) {
      c.beginPath(); c.moveTo(A.X(b.x[0]), A.Y(b.hi[0]));
      for (i = 1; i < b.x.length; i++) c.lineTo(A.X(b.x[i]), A.Y(b.hi[i]));
      for (i = b.x.length - 1; i >= 0; i--) c.lineTo(A.X(b.x[i]), A.Y(b.lo[i]));
      c.closePath(); c.fillStyle = "rgba(" + K.drug + "," + alpha + ")"; c.fill();
      A.line(b.x, b.mean, "rgba(" + K.drug + "," + lineAlpha + ")", 2);
    }
    var saw = st.mode === "patch";
    function lane(alpha, lineAlpha) {
      if (saw) band(alpha, lineAlpha);
      else A.line(b.x, b.hi, "rgba(" + K.drug + "," + lineAlpha + ")", 2);
    }
    lane(0.1, 0.3);
    c.save(); c.beginPath(); c.rect(A.x0, A.y1 - 2, A.X(play) - A.x0, A.y0 - A.y1 + 4); c.clip(); lane(0.3, 1); c.restore();
    var ds = S.RATES.dutyStar;
    A.hline(ds, K.thr, [5, 4], 1.3);
    A.text((saw ? "average needed to close · " : "needed to close · ") + Math.round(ds * 100) + " %",
      A.x1 - 4, A.Y(ds) - 8, K.thr, "right", lbl);

    /* lane 2: front position */
    var rfMax = 0; for (i = 0; i < st.rf.length; i++) if (st.rf[i] > rfMax) rfMax = st.rf[i];
    var ymax = Math.max(st.Rw * 10, rfMax * 10) * 1.12;
    var Bx = Plot.Axes(s, { l: 50, t: T + laneH + gap, r: 14, b: Bm }, [0, xmax], [0, ymax], {
      xlabel: "time (days)", ylabel: "front (mm)", xticks: xt,
      yticks: [[0, "0"], [st.Rw * 10, (st.Rw * 10).toFixed(1)]]
    });
    Bx.frame();
    Bx.hline(st.Rw * 10, K.edge, [4, 4], 1);
    if (rfMax > st.Rw * 1.05) {
      c.fillStyle = "rgba(194,72,79,.06)"; c.fillRect(Bx.x0, Bx.y1, Bx.x1 - Bx.x0, Bx.Y(st.Rw * 10) - Bx.y1);
      Bx.text("above the dashed line: tissue around the wound has also fallen below closure density",
        Bx.x1 - 6, Bx.Y(st.Rw * 10) - 10, "#9C3B42", "right", "500 10.5px Inter, system-ui, sans-serif");
    }
    var ts = [], ys = [];
    for (i = 0; i < st.rf.length; i++) { ts.push(i * FE); ys.push(st.rf[i] * 10); }
    if (st.res.closed) { ts.push(st.endDay); ys.push(0); }
    Bx.line(ts, ys, "#CBD3DA", 2);
    c.save(); c.beginPath(); c.rect(Bx.x0, Bx.y1 - 4, Bx.X(play) - Bx.x0, Bx.y0 - Bx.y1 + 8); c.clip();
    Bx.line(ts, ys, K.front, 3); c.restore();

    // closure marker
    if (st.res.closed) {
      var reached = st.day >= st.endDay - 1e-6, xc = Bx.X(st.endDay);
      c.strokeStyle = reached ? "#1B6B45" : "#B6C0B9"; c.lineWidth = 1.2; c.setLineDash([3, 3]);
      c.beginPath(); c.moveTo(xc, Bx.y0); c.lineTo(xc, Bx.y1); c.stroke(); c.setLineDash([]);
      Bx.text("closed · day " + st.endDay.toFixed(0), xc - 6, Bx.Y(st.Rw * 10) + 11, reached ? "#1B6B45" : C.mute, "right", lbl);
    }

    // playhead through both lanes
    var xp = Bx.X(play);
    c.strokeStyle = "rgba(22,32,27,.3)"; c.lineWidth = 1; c.beginPath(); c.moveTo(xp, A.y1); c.lineTo(xp, Bx.y0); c.stroke();
    Bx.dot(play, cur.rf * 10, K.front, 5);
    A.dot(play, sAt(play), "rgb(" + K.drug + ")", 4);
  }

  function status(cur) {
    var k = cur.k, rf = cur.rf * 10, closed = st.res.closed, done = st.day >= st.endDay - 1e-6, msg;
    var minRf = Math.min.apply(null, st.rf.slice(0, k + 1));
    var duty = st.mode === "patch" ? coverage() / st.iv : 0;
    var centre = Math.round(cur.n[0] * 100), thrPct = Math.round(S.RATES.thr * 100), need = Math.round(S.RATES.dutyStar * 100);
    if (closed && done) {
      msg = "<b>Closed on day " + st.endDay.toFixed(1) + ".</b> The front reached the centre" +
        (st.mode === "patch" ? " after " + patchesBy(st.endDay) + " patches, one every " + st.iv + " h." : ".");
    } else if (!closed && done) {
      if (st.mode === "patch") {
        msg = "<b>Not closed.</b> A patch every " + st.iv + " h gives drug " + Math.round(duty * 100) +
          " % of the time; the wound needs ~" + need + " %. Cells fill in but level off at " + centre +
          " % density, under the " + thrPct + " % that counts as closed.";
      } else {
        msg = st.rf[k] > minRf * 1.25 + 1e-3
          ? "<b>Relapse.</b> The front advanced to " + (minRf * 10).toFixed(1) + " mm, then receded to " + rf.toFixed(1) + " mm once the drug effect faded."
          : "<b>Stalled</b> at " + rf.toFixed(1) + " mm. The wound has not closed by day " + DAYMAX + ".";
      }
    } else {
      msg = "Day " + st.day.toFixed(1) + ": front at <b>" + rf.toFixed(2) + " mm</b>";
      if (st.mode === "patch") msg += " · patch #" + patchesBy(st.day) + " · centre density " +
        centre + " % (closed at " + thrPct + " %)";
      else if (st.mode === "pulse" && st.day > st.T) msg += " (treatment ended, drug effect relaxing)";
    }
    el("i2-status").innerHTML = msg;
  }

  /* ---------- playback ---------- */
  function play(on) {
    st.playing = on; el("i2-play").innerHTML = on ? "&#10074;&#10074; Pause" : "&#9654; Play";
    if (on) {
      if (st.day >= st.endDay - 1e-6) st.day = 0;
      st.parked = false; st.lastPatch = patchesBy(st.day); st.last = performance.now(); requestAnimationFrame(tick);
    }
  }
  // The loop only runs while the widget is on screen.
  function tick(now) {
    if (!st.playing) return;
    if (!st.visible) { st.parked = true; return; }
    var v = speedDays();
    st.day = Math.min(st.endDay, st.day + Math.min(0.1, (now - st.last) / 1000) * v); st.last = now;
    if (st.mode === "patch") {
      var p = patchesBy(st.day);
      // only flash when individual patches are slow enough to see (<= 3 per second)
      if (p !== st.lastPatch && v * 24 / st.iv <= 3) st.flashAt = now;
      st.lastPatch = p;
    }
    el("i2-d").value = st.day;
    draw(now);
    if (st.day >= st.endDay - 1e-6) { play(false); return; }
    requestAnimationFrame(tick);
  }

  function refresh(autoplay) {
    el("i2-T-v").textContent = st.T.toFixed(1) + " d";
    el("i2-iv-v").textContent = st.iv + " h";
    el("i2-rw-v").textContent = (st.Rw * 10).toFixed(1) + " mm";
    el("i2-iv-c").style.display = st.mode === "patch" ? "" : "none";
    el("i2-T-c").style.display = st.mode === "pulse" ? "" : "none";
    simulate(); el("i2-d").value = 0; draw();
    if (autoplay && !reduce) play(true); else play(false);
  }
  var refreshSoon = Plot.debounce(function () { refresh(true); }, 140);

  el("i2-T").addEventListener("input", function () { st.T = parseFloat(this.value); el("i2-T-v").textContent = st.T.toFixed(1) + " d"; refreshSoon(); });
  el("i2-iv").addEventListener("input", function () { st.iv = parseFloat(this.value); el("i2-iv-v").textContent = st.iv + " h"; refreshSoon(); });
  el("i2-rw").addEventListener("input", function () { st.Rw = parseFloat(this.value) / 10; el("i2-rw-v").textContent = (st.Rw * 10).toFixed(1) + " mm"; refreshSoon(); });
  function group(map, key, onPick) {
    Object.keys(map).forEach(function (k) {
      var b = el(map[k]); if (!b) return;
      b.addEventListener("click", function () {
        Object.keys(map).forEach(function (m) { var e = el(map[m]); if (e) { e.classList.toggle("on", m === k); e.setAttribute("aria-pressed", m === k); } });
        onPick(k);
      });
    });
  }
  group({ patch: "i2-pat", pulse: "i2-trt", none: "i2-unt" }, "mode", function (m) { st.mode = m; refresh(true); });
  group({ slow: "i2-sp-slow", normal: "i2-sp-norm", fast: "i2-sp-fast" }, "speed", function (v) { st.speed = v; st.last = performance.now(); });
  el("i2-play").addEventListener("click", function () { play(!st.playing); });
  el("i2-restart").addEventListener("click", function () { st.day = 0; el("i2-d").value = 0; draw(); if (!reduce) play(true); });
  el("i2-d").addEventListener("input", function () {
    // the slider's max is rounded down to its 0.25-day step, so the last notch
    // snaps to the exact closure day -- otherwise scrubbing never shows "Closed"
    var v = parseFloat(this.value);
    play(false); st.day = v >= parseFloat(this.max) - 1e-9 ? st.endDay : v; draw();
  });
  window.addEventListener("resize", Plot.debounce(function () { draw(); }, 120));

  // the 90-day solve only runs once the widget is about a screen away
  Plot.whenNear(el("i2"), function () {
    refresh(false);
    if ("IntersectionObserver" in window) {
      // autoplay once on first view, park the loop when scrolled away, resume when back
      new IntersectionObserver(function (es) {
        st.visible = es[0].isIntersecting;
        if (!st.visible) return;
        if (!st.autoplayed && !reduce) { st.autoplayed = true; if (st.day === 0) play(true); return; }
        if (st.playing && st.parked) { st.parked = false; st.last = performance.now(); requestAnimationFrame(tick); }
      }, { threshold: 0.3 }).observe(el("i2"));
    } else st.visible = true;
  });
})();
