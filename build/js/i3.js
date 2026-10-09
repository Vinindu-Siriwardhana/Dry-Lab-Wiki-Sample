/* Interactive 3 -- scaled dose-response (closed form, no solver). */
(function () {
  "use strict";
  var cv = document.getElementById("i3-cv");
  if (!cv) return;
  var C = Plot.COL, S = Sim, refs = [0.1, 1, 10, 100, 1000];

  function num(id) { return parseFloat(document.getElementById(id).value); }
  function nice(v) { return v >= 100 ? v.toFixed(0) : v >= 10 ? v.toFixed(1) : v >= 1 ? v.toFixed(2) : v.toFixed(3); }

  function draw() {
    var psi = Math.pow(10, num("i3-psi")), a = Math.pow(10, num("i3-a"));
    document.getElementById("i3-psi-v").textContent = nice(psi);
    document.getElementById("i3-a-v").textContent = nice(a);
    var supp = S.suppression(psi, a), th = S.thetaLps(psi, a), p50 = 1 + a;
    document.getElementById("i3-supp").textContent = (supp * 100).toFixed(supp < 0.1 ? 1 : 0) + " %";
    document.getElementById("i3-theta").textContent = th.toFixed(th < 0.01 ? 4 : 3);
    document.getElementById("i3-p50").textContent = nice(p50);

    var s = Plot.setup(cv, 0.62), i, xs = [], psiMin = 1e-2, psiMax = 2e4;
    var ax = Plot.Axes(s, { l: 48, t: 14, r: 16, b: 40 }, [psiMin, psiMax], [0, 1.02], {
      xlog: true, xlabel: "scaled dose  ψ", ylabel: "active receptor, relative to untreated",
      xticks: [[0.01, "0.01"], [0.1, "0.1"], [1, "1"], [10, "10"], [100, "100"], [1e3, "10³"], [1e4, "10⁴"]],
      yticks: [[0, "0"], [0.25, "25%"], [0.5, "50%"], [0.75, "75%"], [1, "100%"]]
    });
    ax.frame();
    for (i = 0; i <= 240; i++) xs.push(psiMin * Math.pow(psiMax / psiMin, i / 240));
    ax.clip();
    refs.forEach(function (r) {
      ax.line(xs, xs.map(function (p) { return S.suppression(p, r); }), "#C8D2DA", 1.4);
    });
    ax.hline(0.5, C.axis, [3, 4], 1);
    ax.unclip();
    ax.clip();
    ax.line(xs, xs.map(function (p) { return S.suppression(p, a); }), C.v14, 3);
    ax.vline(p50, C.thresh, [5, 4], 1.4);
    ax.unclip();
    ax.dot(psi, supp, C.fgf, 6);
    ax.dot(p50, 0.5, C.thresh, 4.5);
    ax.text("ψ₅₀ = 1 + a", Math.min(ax.X(p50) + 8, ax.x1 - 86), ax.Y(0.5) - 14, C.thresh);
    ax.text("grey: a = 0.1, 1, 10, 100, 1000", ax.x1 - 4, ax.y1 + 8, C.mute, "right", "500 10.5px Inter, system-ui, sans-serif");
  }

  ["i3-psi", "i3-a"].forEach(function (id) { Plot.slider(id, draw); });
  window.addEventListener("resize", Plot.debounce(draw, 120));
  draw();
})();
