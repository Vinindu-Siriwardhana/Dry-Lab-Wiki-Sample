/* math.js -- render every .eq / .eqi element with KaTeX.
   The element's textContent is raw LaTeX, so no delimiter scanning runs over
   the prose and nothing in the body text can be mangled by accident. */
(function () {
  function render() {
    if (!window.katex) return;
    var nodes = document.querySelectorAll(".eq, .eqi");
    for (var i = 0; i < nodes.length; i++) {
      var el = nodes[i];
      if (el.getAttribute("data-rendered")) continue;
      var src = el.textContent.trim();
      if (!src) continue;
      try {
        katex.render(src, el, {
          displayMode: el.classList.contains("eq"),
          throwOnError: false,
          strict: false,
          trust: false
        });
        el.setAttribute("data-rendered", "1");
      } catch (e) {
        /* leave the LaTeX visible rather than blanking the equation */
        el.setAttribute("data-rendered", "error");
      }
    }
  }
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", render);
  } else {
    render();
  }
  /* KaTeX is loaded with defer, so re-run once everything has settled. */
  window.addEventListener("load", render);
})();
