/* Interactive 2 -- closure animation: radial Fisher-KPP, n(r,t).
   Three modes: repeated patch changes (each patch covers COV h, the effect
   compartment fades in between), one long pulse, or untreated. */
(function () {
  "use strict";
  var disc = document.getElementById("i2-disc");
  if (!disc) return;
  var C = Plot.COL, S = Sim, prof = document.getElementById("i2-prof"), tl = document.getElementById("i2-time");
  var DAYMAX = 90, FE = 0.25, SDT = 60;
  var st = { mode: "patch", T: 14, iv: 12, Rw: 0.2, res: null, rf: [], s: null, endDay: DAYMAX, day: 0,
             playing: false, last: 0, visible: false, parked: false, drawn: -1, autoplayed: false };
  var COV = null;   // hours one patch keeps tissue above target (500 um alginate, defaults of Interactive 1)
  var reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  function el(id) { return document.getElementById(id); }
  function mix(a, b, t) { return Math.round(a + (b - a) * t); }
  function colour(n) {
    n = Math.max(0, Math.min(1, n));
    var A = [244, 220, 203], B = [168, 207, 154], D = [46, 125, 79], u, f, g;
    if (n < 0.5) { u = n / 0.5; f = A; g = B; } else { u = (n - 0.5) / 0.5; f = B; g = D; }
    return "rgb(" + mix(f[0], g[0], u) + "," + mix(f[1], g[1], u) + "," + mix(f[2], g[2], u) + ")";
  }
  function coverage() {
    if (COV === null) COV = S.twoDomain({ CVtotUM: S.P.C_SOL_UM, Lgel: S.P.L_GEL, RV: 1, kprot: S.P.K_PROT }).hoursAbove;
    return COV;
  }

  function simulate() {
    var o = { Rw: st.Rw, nr: 300, tEndDays: DAYMAX, frameEvery: FE };
    if (st.mode === "patch") { st.s = S.scheduleS(st.iv, coverage(), DAYMAX); o.sGrid = st.s; o.sdt = SDT; }
    else if (st.mode === "pulse") { st.s = null; o.sFn = S.pulseS(st.T); }
    else { st.s = null; o.sFn = function () { return 0; }; }
    var res = S.fisher(o), dr = res.Rinf / 300, thr = S.RATES.thr;
    st.res = res;
    st.rf = res.frames.map(function (n) { return S.frontAt(n, thr, res.r, dr, 300); });
    st.endDay = res.closed ? res.tc / S.DAY : (res.frames.length - 1) * FE;
    el("i2-d").max = Math.floor(st.endDay * 4) / 4;
    st.day = 0;
  }

  function frameIdx() { return Math.max(0, Math.min(st.res.frames.length - 1, Math.round(st.day / FE))); }
  function patchesBy(day) { return Math.floor(day * 24 / st.iv) + 1; }

  function draw() {
    if (!st.res) return;
    var k = frameIdx(), n = st.res.frames[k], r = st.res.r, Rinf = st.res.Rinf, rfk = st.rf[k];
    el("i2-d-v").textContent = st.day.toFixed(1) + " d";
    drawDisc(n, r, Rinf, rfk);
    drawProfile(n, r, Rinf, rfk);
    drawTimeline();
    status(k);
  }

  function drawDisc(n, r, Rinf, rfk) {
    var s = Plot.setup(disc, 1), c = s.c, cx = s.w / 2, cy = s.h / 2 - 8, R = s.w / 2 - 16, i;
    // one radial-gradient fill sampled from n(r) -- the previous version filled
    // 300 stacked discs per frame (~10^8 pixels) and stalled the page
    var g = c.createRadialGradient(cx, cy, 0, cx, cy, R), N = r.length - 1, STOPS = 96;
    for (i = 0; i <= STOPS; i++) g.addColorStop(i / STOPS, colour(n[Math.round(i / STOPS * N)]));
    c.beginPath(); c.fillStyle = g; c.arc(cx, cy, R, 0, 6.2832); c.fill();
    ring(c, cx, cy, st.Rw / Rinf * R, "#8C6A57", [4, 4], 1.3);
    if (rfk > 0) ring(c, cx, cy, rfk / Rinf * R, "#B2182B", [], 2.2);
    c.strokeStyle = "#D5DBE0"; c.lineWidth = 1.5; c.beginPath(); c.arc(cx, cy, R + 3, 0, 6.2832); c.stroke();
    c.font = "600 11px Inter, system-ui, sans-serif"; c.textAlign = "center"; c.fillStyle = "#2B3137";
    c.fillText("top-down view · red ring = front · dashed = original edge", cx, s.h - 8);
  }
  function ring(c, cx, cy, rad, col, dash, w) {
    c.beginPath(); c.strokeStyle = col; c.lineWidth = w; c.setLineDash(dash);
    c.arc(cx, cy, rad, 0, 6.2832); c.stroke(); c.setLineDash([]);
  }

  function drawProfile(n, r, Rinf, rfk) {
    var s = Plot.setup(prof, 1), mmax = Rinf * 10, xs = [], ys = [], i, step = Math.max(1, Math.floor(r.length / 160));
    var tk = []; for (var v = 0; v <= mmax; v += mmax > 8 ? 2 : 1) tk.push([v, String(v)]);
    var ax = Plot.Axes(s, { l: 44, t: 14, r: 12, b: 48 }, [0, mmax], [0, 1.05], {
      xlabel: "radius (mm)", ylabel: "cell density n / K", xticks: tk,
      yticks: [[0, "0"], [0.5, "0.5"], [1, "1"]]
    });
    ax.frame();
    for (i = 0; i < r.length; i += step) { xs.push(r[i] * 10); ys.push(n[i]); }
    ax.line(xs, ys, C.good, 3);
    ax.hline(S.RATES.thr, C.thresh, [5, 4], 1.3);
    ax.vline(st.Rw * 10, "#8C6A57", [4, 4], 1.2);
    if (rfk > 0) ax.dot(rfk * 10, S.RATES.thr, C.fgf, 4.5);
    ax.text("front threshold", ax.x1 - 4, ax.Y(S.RATES.thr) - 9, C.thresh, "right", "500 10.5px Inter, system-ui, sans-serif");
  }

  function drawTimeline() {
    // keep a front that drifts past the original edge (relapse) inside the chart
    var rfMax = Math.max.apply(null, st.rf) * 10;
    var s = Plot.setup(tl, 0.34), ts = [], ys = [], i, xt = [],
        ymax = Math.min(st.res.Rinf * 10, Math.max(st.Rw * 10, rfMax)) * 1.12;
    for (var d = 0; d <= DAYMAX; d += 15) xt.push([d, String(d)]);
    var ax = Plot.Axes(s, { l: 50, t: 10, r: 14, b: 38 }, [0, DAYMAX], [0, ymax], {
      xlabel: "time (days)", ylabel: "front (mm)", xticks: xt,
      yticks: [[0, "0"], [st.Rw * 10, (st.Rw * 10).toFixed(0)]]
    });
    ax.frame();
    ax.hline(st.Rw * 10, "#8C6A57", [4, 4], 1);   // original wound edge
    if (st.mode === "pulse") {
      ax.c.fillStyle = "rgba(77,146,33,.10)";
      ax.c.fillRect(ax.X(0), ax.y1, ax.X(Math.min(st.T, DAYMAX)) - ax.X(0), ax.y0 - ax.y1);
      ax.text("treatment", ax.X(0) + 6, ax.y1 + 11, C.good, "left", "600 10.5px Inter, system-ui, sans-serif");
    } else if (st.mode === "patch" && st.s) {
      // drug-effect saw-tooth s(t), drawn in the top band of the chart (0 at mid-height, 1 at the top)
      var sx = [], sy = [], stepS = Math.max(1, Math.round(0.02 * S.DAY / SDT)), lim = Math.min(st.s.length, Math.ceil(st.endDay * S.DAY / SDT));
      for (i = 0; i < lim; i += stepS) { sx.push(i * SDT / S.DAY); sy.push(ymax * (0.5 + 0.5 * st.s[i])); }
      ax.line(sx, sy, "rgba(77,146,33,.55)", 1.2);
      ax.text("drug effect (each tooth = a new patch)", ax.X(0) + 6, ax.y1 + 11, C.good, "left", "600 10.5px Inter, system-ui, sans-serif");
    }
    for (i = 0; i < st.rf.length; i++) { ts.push(i * FE); ys.push(st.rf[i] * 10); }
    if (st.res.closed) { ts.push(st.endDay); ys.push(0); }
    ax.line(ts, ys, "#C9D3DB", 2);
    var cut = ts.filter(function (d) { return d <= st.day + 1e-6; });
    ax.line(cut, ys.slice(0, cut.length), C.fgf, 3);
    ax.dot(Math.min(st.day, st.endDay), st.rf[frameIdx()] * 10, C.fgf, 5);
  }

  function status(k) {
    var rf = st.rf[k] * 10, closed = st.res.closed, done = st.day >= st.endDay - 1e-6, msg;
    var minRf = Math.min.apply(null, st.rf.slice(0, k + 1));
    var duty = st.mode === "patch" ? coverage() / st.iv : 0;
    var centre = Math.round(st.res.frames[k][0] * 100), thrPct = Math.round(S.RATES.thr * 100);
    if (closed && done) {
      msg = "<b>Closed on day " + st.endDay.toFixed(1) + ".</b> The front reached the centre" +
        (st.mode === "patch" ? " after " + patchesBy(st.endDay) + " patches, one every " + st.iv + " h." : ".");
    } else if (!closed && done) {
      if (st.mode === "patch") {
        msg = "<b>Not closed.</b> A patch every " + st.iv + " h gives drug " + Math.round(duty * 100) +
          " % of the time; the wound needs ~33 %. Cells fill in but level off at " + centre +
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

  function play(on) {
    st.playing = on; el("i2-play").innerHTML = on ? "&#10074;&#10074; Pause" : "&#9654; Play";
    if (on) { if (st.day >= st.endDay - 1e-6) st.day = 0; st.parked = false; st.drawn = -1; st.last = performance.now(); requestAnimationFrame(tick); }
  }
  // The loop only runs while the widget is on screen, and only repaints when
  // the stored frame changes (one frame per 0.25 d), so it costs nothing while
  // the reader scrolls past.
  function tick(now) {
    if (!st.playing) return;
    if (!st.visible) { st.parked = true; return; }
    st.day = Math.min(st.endDay, st.day + (now - st.last) / 1000 * 6); st.last = now;
    el("i2-d").value = st.day;
    var k = frameIdx();
    if (k !== st.drawn || st.day >= st.endDay - 1e-6) { draw(); st.drawn = k; }
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
  var MODES = { patch: "i2-pat", pulse: "i2-trt", none: "i2-unt" };
  function pick(mode) {
    st.mode = mode;
    for (var m in MODES) { el(MODES[m]).classList.toggle("on", m === mode); el(MODES[m]).setAttribute("aria-pressed", m === mode); }
    refresh(true);
  }
  el("i2-pat").addEventListener("click", function () { pick("patch"); });
  el("i2-trt").addEventListener("click", function () { pick("pulse"); });
  el("i2-unt").addEventListener("click", function () { pick("none"); });
  el("i2-play").addEventListener("click", function () { play(!st.playing); });
  el("i2-restart").addEventListener("click", function () { st.day = 0; el("i2-d").value = 0; draw(); if (!reduce) play(true); });
  el("i2-d").addEventListener("input", function () { play(false); st.day = parseFloat(this.value); draw(); });
  window.addEventListener("resize", Plot.debounce(draw, 120));

  // the 90-day solve only runs once the widget is about a screen away
  Plot.whenNear(el("i2"), function () {
    refresh(false);
    if ("IntersectionObserver" in window) {
      // track visibility for the whole session: autoplay once on first view,
      // park the loop when scrolled away, resume it when scrolled back
      new IntersectionObserver(function (es) {
        st.visible = es[0].isIntersecting;
        if (!st.visible) return;
        if (!st.autoplayed && !reduce) { st.autoplayed = true; if (st.day === 0) play(true); return; }
        if (st.playing && st.parked) { st.parked = false; st.last = performance.now(); requestAnimationFrame(tick); }
      }, { threshold: 0.3 }).observe(el("i2"));
    } else st.visible = true;
  });
})();
