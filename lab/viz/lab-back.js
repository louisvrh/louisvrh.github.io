/* "Back to the lab" link for the hosted visualisers.
   The pages are otherwise unchanged exports. The link sits in a 40px band added above the
   page's own content, so at the top of the page it covers nothing. Once you scroll down it
   slides away and comes back when you scroll up (or when it gets keyboard focus). */
(function () {
  'use strict';
  if (document.querySelector('.lab-back')) return;
  if (window.self !== window.top) return; // embedded on the lab page itself: no back link needed

  var BAND = 40;
  var css =
    '.lab-back{position:fixed;top:6px;left:8px;z-index:2147483000;display:inline-flex;align-items:center;gap:6px;' +
    'height:28px;padding:0 12px;border-radius:999px;font:500 12px/1 "IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;' +
    'letter-spacing:0;text-decoration:none;color:#1f2937;background:rgba(255,255,255,.88);border:1px solid rgba(0,0,0,.14);' +
    'box-shadow:0 1px 6px rgba(0,0,0,.08);-webkit-backdrop-filter:blur(8px);backdrop-filter:blur(8px);' +
    'transition:transform .2s ease,opacity .2s ease}' +
    '.lab-back:hover{border-color:#0e7490;color:#0e7490}' +
    '.lab-back:focus-visible{outline:2px solid #0e7490;outline-offset:2px}' +
    '.lab-back.lab-back--away{transform:translateY(-48px);opacity:0;pointer-events:none}' +
    '.lab-back:focus-visible.lab-back--away{transform:none;opacity:1;pointer-events:auto}' +
    '@media (prefers-color-scheme: dark){.lab-back{color:#c9d1d9;background:rgba(8,10,15,.82);border-color:rgba(255,255,255,.14)}' +
    '.lab-back:hover{color:#5fd4e6;border-color:#5fd4e6}.lab-back:focus-visible{outline-color:#5fd4e6}}' +
    '@media (prefers-reduced-motion: reduce){.lab-back{transition:none}}' +
    '@media print{.lab-back{display:none}}';

  function init() {
    var style = document.createElement('style');
    style.textContent = css;
    document.head.appendChild(style);

    // Reserve the band at the top so the link never sits on the page's title or controls.
    var pad = parseFloat(getComputedStyle(document.body).paddingTop) || 0;
    document.body.style.paddingTop = (pad + BAND) + 'px';

    var a = document.createElement('a');
    a.className = 'lab-back';
    a.href = '../../';
    a.innerHTML = '<span aria-hidden="true">←</span> Back to the lab';
    document.body.insertBefore(a, document.body.firstChild);

    var lastY = window.scrollY, ticking = false;
    function update() {
      var y = window.scrollY;
      if (y <= BAND) a.classList.remove('lab-back--away');
      else if (y > lastY + 4) a.classList.add('lab-back--away');
      else if (y < lastY - 4) a.classList.remove('lab-back--away');
      lastY = y;
      ticking = false;
    }
    window.addEventListener('scroll', function () {
      if (!ticking) { ticking = true; requestAnimationFrame(update); }
    }, { passive: true });
  }

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
  else init();
})();
