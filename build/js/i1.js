/* Interactive 1 -- design-space explorer (runs the gel-tissue solver live). */
(function () {
  "use strict";
  var cv = document.getElementById("i1-cv-c");
  if (!cv) return;
  var C = Plot.COL, S = Sim, last = null;
  var DEF = { cv: 738, lg: 500, rv: 0, kp: -4.5, rw: 2, iv: 12 };

  function val(id) { return parseFloat(document.getElementById(id).value); }
  function set(id, v) { document.getElementById(id).value = v; }
  function text(id, s) { document.getElementById(id).textContent = s; }

  function inputs() {
    return { CVtotUM: val("i1-cv"), Lgel: val("i1-lg") * 1e-4, RV: Math.pow(10, val("i1-rv")),
             kprot: Math.pow(10, val("i1-kp")) };
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
    var r = S.twoDomain({ CVtotUM: o.CVtotUM, Lgel: o.Lgel, RV: o.RV, kprot: o.kprot });
    last = r;
    drawChart(r);
    text("i1-hours", r.hoursAbove > 0 ? Plot.fmt(r.hoursAbove, r.hoursAbove < 10 ? 2 : 1) + " h" : "0 h");
    text("i1-peak", Math.round(r.peakUM) + " µM");
    var pct = r.solFrac * 100, tile = document.getElementById("i1-sol-tile");
    text("i1-sol", pct.toFixed(0) + " %");
    text("i1-sol-k", r.overCeiling ? "of the solubility ceiling — CANNOT be formulated" : "of the 738 µM solubility ceiling");
    tile.classList.toggle("tile-bad", r.overCeiling);
    // the longest workable interval depends only on the patch, not on the
    // wound, so it is recomputed only when the patch's cover actually changes
    var covKey = r.hoursAbove.toFixed(3);
    if (covKey !== maxKey) {
      maxKey = covKey; maxH = null;
      text("i1-int", "…"); text("i1-int-k", "computing the longest workable change interval");
      document.getElementById("i1-int-tile").classList.remove("tile-bad");
      searchMax(r);
    }
    protocol();
  }

  /* ---- outputs that depend on the treatment protocol (radius, interval) ---- */
  var maxKey = null, maxH = null, cancelSearch = null, tcTimer = null;

  var searchMax = Plot.debounce(function (r) {
    if (cancelSearch) cancelSearch();
    cancelSearch = S.maxIntervalAsync(0.2, r.hoursAbove, function (h) {
      cancelSearch = null;
      if (r !== last) return;
      maxH = h; showMax(h); protocol();
    }, 24, 200);          // 200-point grid: same interval as 300, half the cost
  }, 220);

  function showMax(h) {
    var tile = document.getElementById("i1-int-tile");
    tile.classList.toggle("tile-bad", h <= 0);
    if (h <= 0) {
      text("i1-int", "none");
      text("i1-int-k", "this patch never holds tissue above target, so no change interval closes a wound");
    } else if (h >= 24) {
      text("i1-int", "≥ 24 h");
      text("i1-int-k", "longest change interval that still closes a wound, for every wound size; even a daily change works");
    } else {
      text("i1-int", "≤ " + h.toFixed(1) + " h");
      text("i1-int-k", "longest change interval that still closes a wound, for every wound size; a daily change never does");
    }
  }

  function protocol() {
    var Rw = val("i1-rw") / 10, iv = val("i1-iv");
    text("i1-rw-v", (Rw * 10).toFixed(1) + " mm");
    text("i1-iv-v", iv.toFixed(1) + " h");
    text("i1-tc", "…");
    text("i1-tc-k", "simulating closure of a " + (Rw * 10).toFixed(1) + " mm wound");
    clearTimeout(tcTimer);
    tcTimer = setTimeout(function () { closureTime(Rw, iv); }, 160);
  }

  function closureTime(Rw, iv) {
    if (!last) return;
    var cov = last.hoursAbove, tile = document.getElementById("i1-tc-tile"),
        mm = (Rw * 10).toFixed(1) + " mm wound", every = "changed every " + iv.toFixed(1) + " h";
    var tc = cov > 0 ? S.closureTimeInterval(Rw, iv, cov, 120, 200) : Infinity;
    var bad = !isFinite(tc);
    tile.classList.toggle("tile-bad", bad);
    if (!bad) {
      text("i1-tc", tc.toFixed(0) + " d");
      text("i1-tc-k", "to close a " + mm + " with the patch " + every);
    } else if (cov > 0 && maxH !== null && iv <= maxH) {
      text("i1-tc", "> 120 d");
      text("i1-tc-k", "a " + mm + " with the patch " + every + " is still open after 120 days");
    } else {
      text("i1-tc", "never");
      text("i1-tc-k", "a " + mm + " with the patch " + every + " does not close: the gaps between patches are too long");
    }
  }

  ["i1-cv", "i1-lg", "i1-rv", "i1-kp"].forEach(function (id) { Plot.slider(id, update); });
  ["i1-rw", "i1-iv"].forEach(function (id) { Plot.slider(id, protocol); });
  document.getElementById("i1-reset").addEventListener("click", function () {
    set("i1-cv", DEF.cv); set("i1-lg", DEF.lg); set("i1-rv", DEF.rv); set("i1-kp", DEF.kp);
    set("i1-rw", DEF.rw); set("i1-iv", DEF.iv);
    update();
  });
  window.addEventListener("resize", Plot.debounce(function () { if (last) drawChart(last); }, 120));
  Plot.whenNear(document.getElementById("i1"), update);
})();
