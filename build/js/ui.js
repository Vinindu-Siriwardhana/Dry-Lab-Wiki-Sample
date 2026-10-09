/* ui.js -- page chrome: progress bar, scrollspy, back-to-top, lightbox, reveal, copy buttons. */
(function () {
  "use strict";
  var doc = document, reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
  function $(s, r) { return (r || doc).querySelector(s); }
  function $$(s, r) { return Array.prototype.slice.call((r || doc).querySelectorAll(s)); }

  // progress bar + back to top + nav shadow
  var bar = doc.createElement("div"); bar.className = "progressbar"; doc.body.appendChild(bar);
  var top = doc.createElement("button"); top.className = "totop"; top.type = "button";
  top.setAttribute("aria-label", "Back to top"); top.innerHTML = "&#8593;"; doc.body.appendChild(top);
  top.addEventListener("click", function () { window.scrollTo({ top: 0, behavior: reduce ? "auto" : "smooth" }); });
  var nav = $("nav.toc");
  // read layout once per frame at most, and draw the bar with a transform
  // (compositor only) instead of animating its width
  var ticking = false, maxScroll = 0;
  function measure() { var h = doc.documentElement; maxScroll = h.scrollHeight - h.clientHeight; }
  function paint() {
    ticking = false;
    var y = window.scrollY;
    bar.style.transform = "scaleX(" + (maxScroll > 0 ? Math.min(1, y / maxScroll) : 0) + ")";
    top.classList.toggle("show", y > 700);
    if (nav) nav.classList.toggle("scrolled", y > 80);
  }
  function onScroll() { if (!ticking) { ticking = true; requestAnimationFrame(paint); } }
  window.addEventListener("scroll", onScroll, { passive: true });
  window.addEventListener("resize", function () { measure(); onScroll(); });
  window.addEventListener("load", function () { measure(); onScroll(); });
  if ("ResizeObserver" in window) new ResizeObserver(function () { measure(); onScroll(); }).observe(doc.body);
  measure(); paint();

  // scrollspy
  // Only the nav strip itself is scrolled sideways to keep the active pill in
  // view. Element.scrollIntoView() must NOT be used here: it also scrolls the
  // page, which yanked the reader backwards at every section boundary.
  var links = $$("nav.toc a"), map = {}, strip = $("nav.toc ul"), current = null;
  links.forEach(function (a) { map[a.getAttribute("href").slice(1)] = a; });
  var secs = Object.keys(map).map(function (id) { return doc.getElementById(id); }).filter(Boolean);
  function activate(id) {
    var a = map[id];
    if (!a || a === current) return;
    if (current) current.classList.remove("active");
    a.classList.add("active"); current = a;
    if (strip && strip.scrollWidth > strip.clientWidth + 1) {
      var left = a.offsetLeft - (strip.clientWidth - a.offsetWidth) / 2;
      strip.scrollTo({ left: Math.max(0, left), behavior: reduce ? "auto" : "smooth" });
    }
  }
  if ("IntersectionObserver" in window && secs.length) {
    var spy = new IntersectionObserver(function (es) {
      es.forEach(function (e) { if (e.isIntersecting) activate(e.target.id); });
    }, { rootMargin: "-30% 0px -60% 0px" });
    secs.forEach(function (s) { spy.observe(s); });
  }

  // reveal on scroll
  if ("IntersectionObserver" in window && !reduce) {
    var items = $$(".claim, figure, .widget, .tablewrap, .keystrip > div, .asks li, .corr li, .filecard");
    items.forEach(function (n) { n.classList.add("reveal"); });
    var rv = new IntersectionObserver(function (es) {
      es.forEach(function (e) { if (e.isIntersecting) { e.target.classList.add("in"); rv.unobserve(e.target); } });
    }, { rootMargin: "0px 0px -8% 0px", threshold: 0.05 });
    items.forEach(function (n) { rv.observe(n); });
  }

  // figure lightbox
  var lb = doc.createElement("div"); lb.className = "lightbox"; lb.setAttribute("role", "dialog"); lb.setAttribute("aria-modal", "true");
  lb.innerHTML = '<button class="lb-x" type="button" aria-label="Close">&times;</button><img alt=""><p></p>';
  doc.body.appendChild(lb);
  var lbImg = $("img", lb), lbCap = $("p", lb);
  var lastFocus = null;
  function closeLb() {
    if (!lb.classList.contains("open")) return;
    lb.classList.remove("open"); doc.body.classList.remove("noscroll");
    if (lastFocus) lastFocus.focus({ preventScroll: true });
  }
  $$(".figbox img").forEach(function (img) {
    img.setAttribute("tabindex", "0"); img.setAttribute("role", "button");
    img.setAttribute("aria-label", "Enlarge figure: " + (img.alt || "").slice(0, 70));
    function open() {
      lbImg.src = img.src; lbImg.alt = img.alt;
      var cap = img.closest("figure") && img.closest("figure").querySelector(".fid");
      lbCap.textContent = cap ? cap.textContent.replace(/ /g, " ") + " Click anywhere or press Esc to close." : "";
      lastFocus = img;
      lb.classList.add("open"); doc.body.classList.add("noscroll"); $(".lb-x", lb).focus({ preventScroll: true });
    }
    img.addEventListener("click", open);
    img.addEventListener("keydown", function (e) { if (e.key === "Enter" || e.key === " ") { e.preventDefault(); open(); } });
  });
  lb.addEventListener("click", closeLb);
  doc.addEventListener("keydown", function (e) { if (e.key === "Escape") closeLb(); });

  // code copy buttons + model-map smooth scroll
  $$(".copybtn").forEach(function (b) {
    b.addEventListener("click", function () {
      var code = b.closest("details").querySelector("code");
      var done = function () { b.textContent = "copied"; setTimeout(function () { b.textContent = "copy"; }, 1600); };
      function sel() { var r = doc.createRange(); r.selectNodeContents(code); var s = window.getSelection(); s.removeAllRanges(); s.addRange(r); }
      try { navigator.clipboard.writeText(code.textContent).then(done, sel); } catch (e) { sel(); }
    });
  });
  $$(".modelmap a").forEach(function (a) {
    a.addEventListener("click", function (e) {
      var t = $(a.getAttribute("href")); if (t) { e.preventDefault(); t.scrollIntoView({ behavior: reduce ? "auto" : "smooth", block: "start" }); }
    });
  });
})();
