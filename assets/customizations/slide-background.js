// Applies sidebar background + content-offset classes based on slide type.
// Title slides use `.left-bg-image` (wide sidebar); content slides use
// `.narrow-left-bg-image` (narrow sidebar, with a horizontal rule under the title).
document.addEventListener('DOMContentLoaded', function () {
  function updateSlideStyles() {
    const currentSlide = Reveal.getCurrentSlide();
    if (!currentSlide) return;

    // Reset body-level background/active classes
    document.body.classList.remove('wide-sidebar-active', 'narrow-sidebar-active');
    document.body.classList.remove('wide-sidebar-bg', 'narrow-sidebar-bg');

    // Remove any existing horizontal rules
    document.querySelectorAll('.horizontal-line').forEach(line => line.remove());

    // Clear content-offset + background classes from all slides to avoid inheritance
    document.querySelectorAll('.slides section').forEach(slide => {
      slide.querySelectorAll('h1, h2, h3, p, ul, ol, .columns, [class*="column"]').forEach(el => {
        el.classList.remove('content-with-wide-sidebar', 'content-with-narrow-sidebar');
      });
      slide.classList.remove('wide-sidebar-bg', 'narrow-sidebar-bg');
    });

    if (currentSlide.classList.contains('left-bg-image')) {
      // Title slide -- wide sidebar
      document.body.classList.add('wide-sidebar-bg');
      document.body.classList.add('wide-sidebar-active');
      currentSlide.querySelectorAll('h1, h2, h3, p, ul, ol, .columns, [class*="column"]').forEach(el => {
        el.classList.add('content-with-wide-sidebar');
      });
    } else if (currentSlide.classList.contains('narrow-left-bg-image')) {
      // Content slide -- narrow sidebar
      document.body.classList.add('narrow-sidebar-bg');
      document.body.classList.add('narrow-sidebar-active');
      currentSlide.querySelectorAll('h1, h2, h3, p, ul, ol, .columns, [class*="column"]').forEach(el => {
        el.classList.add('content-with-narrow-sidebar');
      });
      const horizontalLine = document.createElement('div');
      horizontalLine.className = 'horizontal-line';
      currentSlide.appendChild(horizontalLine);
    }
  }

  Reveal.on('slidechanged', updateSlideStyles);
  Reveal.on('ready', updateSlideStyles);

  const originalLayout = Reveal.layout;
  Reveal.layout = function () {
    originalLayout.apply(this, arguments);
    updateSlideStyles();
  };
});
