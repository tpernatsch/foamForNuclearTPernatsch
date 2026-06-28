# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import datetime

# NOTE: temporarily disable Python API doc builds
# sys.path.insert(0, os.path.abspath("../../pythonapi/foamForNuclear"))

current_year = datetime.datetime.now().year


# -- Project information -----------------------------------------------------

project = "foamForNuclear"
copyright = f"2015-{current_year}, foamForNuclear Team"
author = "foamForNuclear Team"


# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.napoleon",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "myst_parser",
    "sphinxmermaid",
]

extensions += ["sphinxcontrib.bibtex"]
bibtex_bibfiles = [
    "usersguide/OFFBEAT/references.bib",
]

intersphinx_mapping = {
    "rtd": ("https://docs.readthedocs.com/platform/stable/", None),
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

myst_enable_extensions = [
    "amsmath",
    "attrs_inline",
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    # "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

templates_path = ["_templates"]

# -- Options for EPUB output
epub_show_urls = "footnote"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# exclude_patterns += [
#     "pythonapi/**",
#     "**/pythonapi/**",
#     "**/pyapi/**",
#     "api/**",
#     "**/api/**",
# ]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_title = "foamForNuclear Documentation"
html_logo = "../logo/ffn-logo-bold-italic.svg"
html_favicon = "../logo/ffn-logo.svg"

html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": True,
}

# Don't ship the raw .rst/.md sources or "view source" links — saves the
# _sources/ tree in the published site.
html_copy_source = False
html_show_sourcelink = False


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]
html_css_files = ["theme_overrides.css"]

source_suffix = {
    ".rst": "restructuredtext",
    ".txt": "markdown",
    ".md": "markdown",
}

# ---------------------------------------------------------------------------
# Dual-build: standard vs AI-enhanced
# ---------------------------------------------------------------------------

from pathlib import Path as _Path

_CONF_DIR = _Path(__file__).resolve().parent        # documentation/sphinx/
_REPO_ROOT = _CONF_DIR.parent.parent                # repo root
_AI_CONTENT_DIR = _REPO_ROOT / "tools" / "docAI" / "ai_content"

AI_ENHANCED = os.getenv("SPHINX_AI", "0") == "1"

if AI_ENHANCED:
    html_title = "foamForNuclear Documentation (AI-enhanced)"
    html_context = {"ai_enhanced": True}
else:
    html_context = {"ai_enhanced": False}

# Note shown on AI-polished usersguide pages so readers can tell which pages
# differ from the standard documentation.
_AI_NOTE_RST = """
.. note::

   The language of this page was automatically polished by an AI model.
   Use the *Standard docs* switcher button to view the original page.
"""

_AI_NOTE_MD = """
```{note}
The language of this page was automatically polished by an AI model.
Use the *Standard docs* switcher button to view the original page.
```
"""


def _inject_ai_note(text, suffix):
    """Append the AI-polished note at the end of the page."""
    note = _AI_NOTE_MD if suffix in (".md", ".txt") else _AI_NOTE_RST
    return text.rstrip() + "\n\n" + note + "\n"


# Inject AI content at read-time for the AI build (no file copies needed).
# For each usersguide page, if an AI-polished version exists, replace the
# source string; otherwise fall through silently to the original.
def _overlay_ai_usersguide(app, docname, source):
    if not docname.startswith("usersguide/"):
        return
    # The parser was chosen from the suffix of the real source file, so the
    # overlay must use the ai_content file with that same suffix — injecting
    # another format would garble the page (e.g. RST parsed as Markdown).
    src_suffix = _Path(app.env.doc2path(docname)).suffix
    ai_path = _AI_CONTENT_DIR / "documentation" / "sphinx" / (docname + src_suffix)
    if ai_path.exists():
        source[0] = _inject_ai_note(
            ai_path.read_text(encoding="utf-8"), src_suffix
        )


def setup(app):
    app.add_js_file('filter-table.js')
    app.add_js_file('version_switch.js')
    if AI_ENHANCED:
        app.connect('source-read', _overlay_ai_usersguide)