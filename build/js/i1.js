/* Interactive 1 -- design-space explorer (runs the gel-tissue solver live). */
(function () {
  "use strict";
  var cv = document.getElementById("i1-cv-c");
  if (!cv) return;
  var C = Plot.COL, S = Sim, last = null;
  var DEF = { cv: 738, lg: 150, rv: 0, kp: -4.5, rw: 2 };

  function val(id) { return parseFloat(document.getElementById(id).value); }
  function set(id, v) { document.getElementById(id).value = v; }
  function text(id, s) { document.getElementById(id).textContent = s; }

  function inputs() {
    return { CVtotUM: val("i1-cv"), Lgel: val("i1-lg") * 1e-4, RV: Math.pow(10, val("i1-rv")),
             kprot: Math.pow(10, val("i1-kp")), Rw: val("i1-rw") / 10 };
  }

  function drawChart(r) {
    var s = Plot.setup(cv, 0.52), i, hrs = r.t.map(function (t) { return t / 3600; });
    var tmax = 72, ymax = Math.max(40, r.peakUM * 1.15), target = S.P.TARGET_UM;
    var ax = Plot.Axes(s, { l: 50, t: 12, r: 14, b: 40 }, [0.01, tmax], [0, ymax], {
      xlog: true, xlabel: "time after application", ylabel: "tissue V14, depth-averaged (µM)",
      xticks: [[0.0167, "1 min"], [0.1667, "10 min"], [1, "1 h"], [6, "6 h"], [24, "1 d"], [72, "3 d"]],
      yticks: niceTicks(ymax)
    });
    ax.frame();
    // shade the time above target
    ax.c.save(); ax.c.beginPath(); ax.c.rect(ax.x0, ax.y1, ax.x1 - ax.x0, ax.Y(target) - ax.y1); ax.c.clip();
    ax.c.beginPath(); ax.c.fillStyle = "rgba(33,102,172,.20)";
    var started = false;
    for (i = 0; i < hrs.length; i++) {
      if (hrs[i] < 0.01) continue;
      var px = ax.X(hrs[i]), py = ax.Y(r.avg[i] / S.UM);
      if (!started) { ax.c.moveTo(px, ax.Y(0)); ax.c.lineTo(px, py); started = true; } else ax.c.lineTo(px, py);
    }
    ax.c.lineTo(ax.X(hrs[hrs.length - 1]), ax.Y(0)); ax.c.closePath(); ax.c.fill(); ax.c.restore();
    ax.clip();
    ax.line(hrs, r.avg.map(function (v) { return v / S.UM; }), C.v14, 3);
    ax.unclip();
    ax.hline(target, C.thresh, [6, 4], 1.6);
    ax.text("25 µM target", ax.x1 - 6, ax.Y(target) - 10, C.thresh, "right");
    if (r.peakUM > 0) ax.dot(r.peakT / 3600, r.peakUM, C.fgf, 5);
  }

  function sci(v) {
    var e = Math.floor(Math.log10(v)), m = v / Math.pow(10, e), sup = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
    if (m >= 9.95) { m = 1; e += 1; }
    return m.toFixed(1) + "×10" + String(e).split("").map(function (ch) { return sup[ch]; }).join("");
  }

  function niceTicks(ymax) {
    var step = ymax > 400 ? 100 : ymax > 200 ? 50 : ymax > 100 ? 25 : ymax > 50 ? 10 : 5, t = [];
    for (var v = 0; v <= ymax; v += step) t.push([v, String(v)]);
    return t;
  }

  function update() {
    var o = inputs();
    text("i1-cv-v", Math.round(o.CVtotUM) + " µM");
    text("i1-lg-v", Math.round(o.Lgel * 1e4) + " µm");
    text("i1-rv-v", o.RV < 10 ? o.RV.toFixed(1) : o.RV.toFixed(0));
    text("i1-kp-v", sci(o.kprot) + " s⁻¹");
    text("i1-rw-v", (o.Rw * 10).toFixed(1) + " mm");
    var r = S.twoDomain({ CVtotUM: o.CVtotUM, Lgel: o.Lgel, RV: o.RV, kprot: o.kprot });
    last = r;
    drawChart(r);
    text("i1-hours", r.hoursAbove > 0 ? Plot.fmt(r.hoursAbove, r.hoursAbove < 10 ? 2 : 1) + " h" : "0 h");
    text("i1-peak", Math.round(r.peakUM) + " µM");
    var pct = r.solFrac * 100, tile = document.getElementById("i1-sol-tile");
    text("i1-sol", pct.toFixed(0) + " %");
    text("i1-sol-k", r.overCeiling ? "of the solubility ceiling — CANNOT be formulated" : "of the 738 µM solubility ceiling");
    tile.classList.toggle("tile-bad", r.overCeiling);
    text("i1-int", "…"); text("i1-int-k", "computing the patch-change interval");
    document.getElementById("i1-int-tile").classList.remove("tile-bad");
    closure(o, r);
  }

  var closure = Plot.debounce(function (o, r) {
    if (r !== last) return;
    var h = S.maxInterval(o.Rw, r.hoursAbove, 24, 300), tile = document.getElementById("i1-int-tile");
    var mm = (o.Rw * 10).toFixed(1) + " mm wound";
    if (h <= 0) {
      text("i1-int", "never");
      text("i1-int-k", "this patch does not reach the target, so no change interval closes a " + mm);
      tile.classList.add("tile-bad");
    } else if (h >= 24) {
      text("i1-int", "≥ 24 h");
      text("i1-int-k", "longest interval that closes a " + mm + "; even a daily change works");
    } else {
      text("i1-int", "≤ " + h.toFixed(1) + " h");
      text("i1-int-k", "longest patch-change interval that closes a " + mm + "; a daily change does not close it");
    }
  }, 220);

  ["i1-cv", "i1-lg", "i1-rv", "i1-kp", "i1-rw"].forEach(function (id) { Plot.slider(id, update); });
  document.getElementById("i1-reset").addEventListener("click", function () {
    set("i1-cv", DEF.cv); set("i1-lg", DEF.lg); set("i1-rv", DEF.rv); set("i1-kp", DEF.kp); set("i1-rw", DEF.rw);
    update();
  });
  window.addEventListener("resize", Plot.debounce(function () { if (last) drawChart(last); }, 120));
  update();
})();
