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
    # NOTE: temporarily disable Python API doc builds
    # "sphinx.ext.autodoc",
    # "sphinx.ext.autosummary",
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

# NOTE: temporarily exclude Python API pages from the build
exclude_patterns += [
    "pythonapi/**",
    "**/pythonapi/**",
    "**/pyapi/**",
    "api/**",
    "**/api/**",
]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_rtd_theme"
html_title = "foamForNuclear Documentation"
html_logo = "../logo/ffn-logo-bold-italic.svg"
html_favicon = "../logo/ffn-logo.svg"

html_theme_options = {
    "navigation_depth": -1,
    "collapse_navigation": False,
}


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

def setup(app):
    app.add_js_file('filter-table.js')