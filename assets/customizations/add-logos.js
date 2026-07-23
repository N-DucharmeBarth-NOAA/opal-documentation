// Injects up to two logos into the top-right of every slide.
// Reads paths from window.opalLogos (set in the deck header). logo2 is optional.
window.addEventListener('load', function () {
  const logos = window.opalLogos || {};
  const logoDiv = document.createElement('div');
  logoDiv.className = 'multiple-logos';

  function addLogo(src, alt) {
    if (!src) return;
    const img = document.createElement('img');
    img.src = src;
    img.alt = alt || '';
    logoDiv.appendChild(img);
  }

  addLogo(logos.logo1, 'opal logo');
  addLogo(logos.logo2, 'partner logo');

  document.querySelector('.reveal').appendChild(logoDiv);

  function updateLogoPosition() {
    logoDiv.style.display = 'flex';
  }

  Reveal.on('slidechanged', updateLogoPosition);
  Reveal.on('ready', updateLogoPosition);
});
