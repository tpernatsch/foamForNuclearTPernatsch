(function () {
    "use strict";

    // FFN_DOC_VARIANT is injected by the Sphinx layout template ("standard" or "ai").
    // Fall back to URL-based detection for locally-served builds without the variable.
    function getVariant() {
        if (window.FFN_DOC_VARIANT) return window.FFN_DOC_VARIANT;
        var m = window.location.pathname.match(/\/(standard|ai)(\/|$)/);
        return m ? m[1] : null;
    }

    function switchVariant() {
        var current = getVariant();
        if (!current) return;
        var target = current === "standard" ? "ai" : "standard";
        // Rewrite URL prefix if present, otherwise just swap the variant segment
        var href = window.location.href;
        if (href.indexOf("/" + current + "/") !== -1) {
            window.location.href = href.replace("/" + current + "/", "/" + target + "/");
        } else {
            // Served at root — no path prefix to rewrite; link goes to sibling dir
            window.location.href = "../" + target + "/";
        }
    }

    document.addEventListener("DOMContentLoaded", function () {
        var current = getVariant();
        if (!current) return;

        var isAI = current === "ai";
        var btn = document.createElement("a");
        btn.id = "version-switch-btn";
        btn.href = "#";
        btn.title = isAI
            ? "Switch to the standard documentation"
            : "Switch to the AI-enhanced documentation";
        btn.innerHTML = isAI
            ? "<span class='vs-icon'>&#128196;</span> Standard docs"
            : "<span class='vs-icon'>&#10024;</span> AI-enhanced docs";
        btn.addEventListener("click", function (e) {
            e.preventDefault();
            switchVariant();
        });

        var target = document.querySelector(".wy-side-nav-search");
        if (target) {
            target.appendChild(btn);
        }
    });
})();
