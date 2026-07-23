// Applies the opal brand font and colour to standalone occurrences of "opal".
// Website pages currently run the equivalent logic from open-in-new-tab.html;
// this standalone file is loaded by RevealJS presentations.
document.addEventListener("DOMContentLoaded", function () {
  const root = document.querySelector(".reveal") || document.body;
  const skipTags = new Set(["SCRIPT", "STYLE", "NOSCRIPT", "TEXTAREA", "PRE"]);
  const opalRegex = /\b(opal)\b/gi;

  const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
    acceptNode(node) {
      if (!node.nodeValue || !node.nodeValue.toLowerCase().includes("opal")) {
        return NodeFilter.FILTER_REJECT;
      }

      const parent = node.parentElement;
      if (!parent || skipTags.has(parent.tagName)) {
        return NodeFilter.FILTER_REJECT;
      }

      if (parent.closest(".opal-brand")) {
        return NodeFilter.FILTER_REJECT;
      }

      if (parent.namespaceURI !== "http://www.w3.org/1999/xhtml") {
        return NodeFilter.FILTER_REJECT;
      }

      return NodeFilter.FILTER_ACCEPT;
    }
  });

  const textNodes = [];
  while (walker.nextNode()) {
    textNodes.push(walker.currentNode);
  }

  textNodes.forEach(function (node) {
    const text = node.nodeValue;
    const fragment = document.createDocumentFragment();
    let lastIndex = 0;
    let match;
    opalRegex.lastIndex = 0;

    while ((match = opalRegex.exec(text)) !== null) {
      if (match.index > lastIndex) {
        fragment.appendChild(document.createTextNode(text.slice(lastIndex, match.index)));
      }

      const span = document.createElement("span");
      span.className = "opal-brand";
      span.textContent = match[0];
      fragment.appendChild(span);
      lastIndex = opalRegex.lastIndex;
    }

    if (lastIndex < text.length) {
      fragment.appendChild(document.createTextNode(text.slice(lastIndex)));
    }

    node.parentNode.replaceChild(fragment, node);
  });
});
