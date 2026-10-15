/* Interactive 3 -- inflammation explorer: Model 3 occupancy -> Model 4 p65 -> Model 5 nitrite.
   Closed forms plus one explicit ODE; constants mirror code/params.py and New Graphs/model. */
(function () {
  "use strict";
  var cv = document.getElementById("i3-c");
  if (!cv) return;
  var C = Plot.COL, S = Sim;
  var KD_L = 0.065, GAMMA = 0.1, H = 2.0;                       // uM; register gamma, h
  var NO = { N: 1e6, q: 5e-19, K: 0.3, p: 2, kturn: Math.LN2 / (4 * 3600), phi: 0.9, p65lps: 0.85 };
  var DEF = { cv: Math.log10(145), l: 0, kd: Math.log10(10.08) };

  function el(id) { return document.getElementById(id); }
  function val(id) { return Math.pow(10, parseFloat(el(id).value)); }
  function nice(v) { return v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2) : v.toFixed(3); }

  /* accumulated nitrite (uM) after 24 h for a normalised p65 level (Model 5) */
  function nitrite24(p65) {
    var drive = NO.N * NO.q * Math.pow(p65, NO.p) / (Math.pow(NO.K, NO.p) + Math.pow(p65, NO.p));
    var dt = 60, R = 0, Rn, no2 = 0;
    for (var t = 0; t < 86400; t += dt) {
      Rn = R + dt * NO.kturn * (drive - R);
      no2 += dt * NO.phi * 0.5 * (R + Rn);
      R = Rn;
    }
    return no2 / S.UM;
  }
  function chain(cvUM, a, kd) {
    var psi = cvUM / kd, pi0 = S.p65(0, a, GAMMA, H), pr = S.p65(psi, a, GAMMA, H) / pi0;
    return {
      lps: S.thetaLps(psi, a), supp: S.suppression(psi, a), p65: pr,
      no2: nitrite24(NO.p65lps * pr)
    };
  }

  function draw() {
    var cvUM = val("i3-cv"), L = val("i3-l"), kd = val("i3-kd"), a = L / KD_L;
    el("i3-cv-v").textContent = nice(cvUM) + " µM";
    el("i3-l-v").textContent = nice(L) + " µM  (a = " + nice(a) + ")";
    el("i3-kd-v").textContent = nice(kd) + " µM";
    var r = chain(cvUM, a, kd), r0 = chain(0, a, kd);
    el("i3-lps").textContent = (r.lps * 100).toFixed(r.lps < 0.1 ? 1 : 0) + " %";
    el("i3-supp").textContent = (r.supp * 100).toFixed(r.supp < 0.1 ? 1 : 0) + " %";
    el("i3-p65").textContent = (r.p65 * 100).toFixed(r.p65 < 0.1 ? 1 : 0) + " %";
    el("i3-no2").textContent = r.no2.toFixed(r.no2 < 10 ? 1 : 0) + " µM";
    el("i3-no2-k").textContent = "predicted nitrite after 24 h (LPS alone: " + r0.no2.toFixed(0) + " µM)";

    var s = Plot.setup(cv, 0.62), xs = [], i, lo = 0.1, hi = 5000;
    var ax = Plot.Axes(s, { l: 48, t: 14, r: 16, b: 40 }, [lo, hi], [0, 1.02], {
      xlog: true, xlabel: "V14 concentration (µM)", ylabel: "relative to LPS alone",
      xticks: [[0.1, "0.1"], [1, "1"], [10, "10"], [100, "100"], [1000, "1000"]],
      yticks: [[0, "0"], [0.25, "25%"], [0.5, "50%"], [0.75, "75%"], [1, "100%"]]
    });
    ax.frame();
    for (i = 0; i <= 160; i++) xs.push(lo * Math.pow(hi / lo, i / 160));
    var rows = xs.map(function (x) { return chain(x, a, kd); });
    ax.clip();
    [[25, C.ink], [145, C.thresh], [738, C.mute]].forEach(function (m) { ax.vline(m[0], m[1], [3, 4], 1.1); });
    ax.line(xs, rows.map(function (q) { return q.supp; }), C.v14, 3);
    ax.line(xs, rows.map(function (q) { return q.p65; }), C.bad, 2, [7, 5]);
    ax.line(xs, rows.map(function (q) { return q.no2 / r0.no2; }), C.thresh, 2.4);
    ax.unclip();
    ax.dot(cvUM, r.supp, C.v14, 5.5);
    ax.dot(cvUM, r.p65, C.bad, 5);
    ax.dot(cvUM, r.no2 / r0.no2, C.thresh, 5);
    var f = "600 11px Inter, system-ui, sans-serif";
    ax.text("active receptor", ax.x0 + 8, ax.y0 - 42, C.v14, "left", f);
    ax.text("nuclear p65 (dashed)", ax.x0 + 8, ax.y0 - 26, C.bad, "left", f);
    ax.text("24 h nitrite", ax.x0 + 8, ax.y0 - 10, C.thresh, "left", f);
  }

  ["i3-cv", "i3-l", "i3-kd"].forEach(function (id) { Plot.slider(id, draw); });
  el("i3-reset").addEventListener("click", function () {
    el("i3-cv").value = DEF.cv; el("i3-l").value = DEF.l; el("i3-kd").value = DEF.kd; draw();
  });
  window.addEventListener("resize", Plot.debounce(draw, 120));
  draw();
})();
