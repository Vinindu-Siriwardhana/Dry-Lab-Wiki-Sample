/* plot.js -- tiny canvas plotting helpers shared by the interactives. */
var Plot = (function () {
  "use strict";
  var COL = {
    v14: "#2166AC", v14l: "#92C5DE", fgf: "#B2182B", thresh: "#E08214", good: "#4D9221",
    bad: "#762A83", ink: "#2B3137", mute: "#6B747D", grid: "#E6EAED", axis: "#9AA3AB"
  };

  /* Size the canvas from its CSS width. The backing store and the CSS height are
     only written when they actually change: rewriting them on every animation
     frame forced a full page re-layout 60 times a second and made scrolling
     stutter while the closure animation played. */
  function setup(canvas, aspect) {
    var dpr = window.devicePixelRatio || 1, w = Math.max(240, Math.round(canvas.clientWidth)),
        h = Math.round(w * aspect), bw = Math.round(w * dpr), bh = Math.round(h * dpr);
    if (canvas._h !== h) { canvas.style.height = h + "px"; canvas._h = h; }
    if (canvas.width !== bw || canvas.height !== bh) { canvas.width = bw; canvas.height = bh; }
    var c = canvas.getContext("2d");
    c.setTransform(dpr, 0, 0, dpr, 0, 0);
    c.clearRect(0, 0, w, h);
    return { c: c, w: w, h: h };
  }

  /* box = {l,t,r,b} pixel margins; xr,yr = [min,max]; o.xlog / o.ylog */
  function Axes(s, box, xr, yr, o) {
    o = o || {};
    var a = { c: s.c, o: o, x0: box.l, x1: s.w - box.r, y0: s.h - box.b, y1: box.t };
    function tf(v, r, log) { return log ? Math.log(v / r[0]) / Math.log(r[1] / r[0]) : (v - r[0]) / (r[1] - r[0]); }
    a.X = function (v) { return a.x0 + tf(v, xr, o.xlog) * (a.x1 - a.x0); };
    a.Y = function (v) { return a.y0 + tf(v, yr, o.ylog) * (a.y1 - a.y0); };
    a.invX = function (px) {
      var f = (px - a.x0) / (a.x1 - a.x0);
      return o.xlog ? xr[0] * Math.pow(xr[1] / xr[0], f) : xr[0] + f * (xr[1] - xr[0]);
    };
    a.frame = function () {
      var c = a.c, i, t;
      c.font = "11px Inter, system-ui, sans-serif"; c.lineWidth = 1;
      c.strokeStyle = COL.grid; c.fillStyle = COL.mute;
      c.textAlign = "center"; c.textBaseline = "top";
      for (i = 0; i < (o.xticks || []).length; i++) {
        t = o.xticks[i]; var x = a.X(t[0]);
        c.beginPath(); c.moveTo(x, a.y0); c.lineTo(x, a.y1); c.stroke();
        c.fillText(t[1], x, a.y0 + 6);
      }
      c.textAlign = "right"; c.textBaseline = "middle";
      for (i = 0; i < (o.yticks || []).length; i++) {
        t = o.yticks[i]; var y = a.Y(t[0]);
        c.beginPath(); c.moveTo(a.x0, y); c.lineTo(a.x1, y); c.stroke();
        c.fillText(t[1], a.x0 - 7, y);
      }
      c.strokeStyle = COL.axis;
      c.beginPath(); c.moveTo(a.x0, a.y1); c.lineTo(a.x0, a.y0); c.lineTo(a.x1, a.y0); c.stroke();
      c.fillStyle = COL.ink; c.font = "600 11.5px Inter, system-ui, sans-serif";
      if (o.xlabel) { c.textAlign = "center"; c.textBaseline = "bottom"; c.fillText(o.xlabel, (a.x0 + a.x1) / 2, s.h - 2); }
      if (o.ylabel) {
        c.save(); c.translate(13, (a.y0 + a.y1) / 2); c.rotate(-Math.PI / 2);
        c.textAlign = "center"; c.textBaseline = "middle"; c.fillText(o.ylabel, 0, 0); c.restore();
      }
    };
    a.clip = function () { a.c.save(); a.c.beginPath(); a.c.rect(a.x0, a.y1, a.x1 - a.x0, a.y0 - a.y1); a.c.clip(); };
    a.unclip = function () { a.c.restore(); };
    a.line = function (xs, ys, color, w, dash) {
      var c = a.c, started = false;
      c.beginPath(); c.strokeStyle = color; c.lineWidth = w || 2; c.setLineDash(dash || []);
      for (var i = 0; i < xs.length; i++) {
        if (!(xs[i] > 0) && o.xlog) continue;
        var px = a.X(xs[i]), py = a.Y(ys[i]);
        if (!started) { c.moveTo(px, py); started = true; } else c.lineTo(px, py);
      }
      c.stroke(); c.setLineDash([]);
    };
    a.hline = function (v, color, dash, w) {
      var y = a.Y(v), c = a.c; c.beginPath(); c.strokeStyle = color; c.lineWidth = w || 1.5;
      c.setLineDash(dash || [5, 4]); c.moveTo(a.x0, y); c.lineTo(a.x1, y); c.stroke(); c.setLineDash([]);
    };
    a.vline = function (v, color, dash, w) {
      var x = a.X(v), c = a.c; c.beginPath(); c.strokeStyle = color; c.lineWidth = w || 1.5;
      c.setLineDash(dash || [5, 4]); c.moveTo(x, a.y0); c.lineTo(x, a.y1); c.stroke(); c.setLineDash([]);
    };
    a.dot = function (x, y, color, r) {
      var c = a.c; c.beginPath(); c.fillStyle = "#fff"; c.arc(a.X(x), a.Y(y), (r || 5) + 2, 0, 6.2832); c.fill();
      c.beginPath(); c.fillStyle = color; c.arc(a.X(x), a.Y(y), r || 5, 0, 6.2832); c.fill();
    };
    a.text = function (str, x, y, color, align, font) {
      var c = a.c; c.fillStyle = color || COL.ink; c.font = font || "600 11.5px Inter, system-ui, sans-serif";
      c.textAlign = align || "left"; c.textBaseline = "middle"; c.fillText(str, x, y);
    };
    return a;
  }

  /* bind a <input type=range> to a value formatter; fn(value) is called on input */
  function slider(id, onInput) {
    var el = document.getElementById(id);
    if (!el) return null;
    el.addEventListener("input", function () { onInput(parseFloat(el.value)); });
    return el;
  }

  function fmt(v, d) {
    if (!isFinite(v)) return "∞";
    return v.toFixed(d == null ? 1 : d);
  }

  function debounce(fn, ms) {
    var h = null;
    return function () { var a = arguments; clearTimeout(h); h = setTimeout(function () { fn.apply(null, a); }, ms); };
  }

  /* run fn once, when `el` comes within ~1 screen of the viewport, and in an
     idle moment, so widget start-up never blocks the first scroll */
  function whenNear(el, fn) {
    var go = function () {
      var ric = window.requestIdleCallback || function (f) { return setTimeout(f, 1); };
      ric(fn, { timeout: 400 });
    };
    if (!el || !("IntersectionObserver" in window)) { go(); return; }
    var io = new IntersectionObserver(function (es) {
      if (es[0].isIntersecting) { io.disconnect(); go(); }
    }, { rootMargin: "900px 0px" });
    io.observe(el);
  }

  return { COL: COL, setup: setup, Axes: Axes, slider: slider, fmt: fmt, debounce: debounce, whenNear: whenNear };
})();
