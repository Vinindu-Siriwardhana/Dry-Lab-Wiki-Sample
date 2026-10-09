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
  function onScroll() {
    var h = doc.documentElement, max = h.scrollHeight - h.clientHeight, y = window.scrollY || h.scrollTop;
    bar.style.width = (max > 0 ? Math.min(100, y / max * 100) : 0) + "%";
    top.classList.toggle("show", y > 700);
    if (nav) nav.classList.toggle("scrolled", y > 80);
  }
  window.addEventListener("scroll", onScroll, { passive: true }); onScroll();

  // scrollspy
  var links = $$("nav.toc a"), map = {};
  links.forEach(function (a) { map[a.getAttribute("href").slice(1)] = a; });
  var secs = Object.keys(map).map(function (id) { return doc.getElementById(id); }).filter(Boolean);
  if ("IntersectionObserver" in window && secs.length) {
    var spy = new IntersectionObserver(function (es) {
      es.forEach(function (e) {
        if (e.isIntersecting) {
          links.forEach(function (l) { l.classList.remove("active"); });
          var a = map[e.target.id]; if (a) { a.classList.add("active"); a.scrollIntoView({ block: "nearest", inline: "center" }); }
        }
      });
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
  function closeLb() { lb.classList.remove("open"); doc.body.classList.remove("noscroll"); }
  $$(".figbox img").forEach(function (img) {
    img.setAttribute("tabindex", "0"); img.setAttribute("role", "button");
    img.setAttribute("aria-label", "Enlarge figure: " + (img.alt || "").slice(0, 70));
    function open() {
      lbImg.src = img.src; lbImg.alt = img.alt;
      var cap = img.closest("figure") && img.closest("figure").querySelector(".fid");
      lbCap.textContent = cap ? cap.textContent.replace(/ /g, " ") + " Click anywhere or press Esc to close." : "";
      lb.classList.add("open"); doc.body.classList.add("noscroll"); $(".lb-x", lb).focus();
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
