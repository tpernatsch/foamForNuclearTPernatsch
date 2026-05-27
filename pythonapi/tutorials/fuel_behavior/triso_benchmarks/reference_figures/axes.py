"""
Axis scales and plot style for the TECDOC benchmark comparison figures.

AXES
    Keyed by SVG stem (filename without extension).
    xlim in units of 10²⁵ n/m², ylim in MPa.
    Needed to denormalise the [0, 1] coordinates from extractSVGPlotData.py.

SERIES_STYLE
    Maps legend label (as it appears in the SVG / XLSX sheet name) to
    matplotlib kwargs.  Used by every post.py comparison plot so that the
    same code always gets the same colour.  Unrecognised labels fall back
    to the default matplotlib cycle.

OFFBEAT_STYLE
    Kwargs for the OFFBEAT curve, always plotted last in red.
"""

AXES = {
    # ── cases 4a ────────────────────────────────────────────────────────────
    "case_4a_IPyC_tangential_stress":             dict(xlim=(0, 3.0), ylim=( -100, 1000), ytick_step=100),
    "case_4a_SiC_tangential_stress":              dict(xlim=(0, 3.0), ylim=(-1000,  200), ytick_step=200),
    "case_4a_interface_radial_stress": dict(xlim=(0, 3.0), ylim=(  -40,  180), ytick_step= 20),

    # ── cases 4b ────────────────────────────────────────────────────────────
    "case_4b_IPyC_tangential_stress":             dict(xlim=(0, 3.0), ylim=(  -30,   15), ytick_step=  5),
    "case_4b_SiC_tangential_stress":              dict(xlim=(0, 3.0), ylim=(  -25,  150), ytick_step= 25),
    "case_4b_interface_radial_stress": dict(xlim=(0, 3.0), ylim=(  -30,    0), ytick_step=  5),

    # ── cases 4c ────────────────────────────────────────────────────────────
    "case_4c_IPyC_tangential_stress":             dict(xlim=(0, 3.0), ylim=(   -5,   30), ytick_step=  5),
    "case_4c_SiC_tangential_stress":              dict(xlim=(0, 3.0), ylim=(  -25,  125), ytick_step= 25),
    "case_4c_interface_radial_stress": dict(xlim=(0, 3.0), ylim=(  -20,    0), ytick_step=  5),

    # ── cases 4d ────────────────────────────────────────────────────────────
    "case_4d_IPyC_tangential_stress":             dict(xlim=(0, 3.0), ylim=(    0,  175), ytick_step= 25),
    "case_4d_SiC_tangential_stress":              dict(xlim=(0, 3.0), ylim=(  -60,  120), ytick_step= 20),
    "case_4d_interface_radial_stress": dict(xlim=(0, 3.0), ylim=(  -20,   10), ytick_step=  5),

    # ── cases 5–7 ───────────────────────────────────────────────────────────
    "case_5_IPyC_tangential_stress": dict(xlim=(0, 3.0), ylim=( -50,  250), ytick_step= 50),
    "case_5_SiC_tangential_stress":  dict(xlim=(0, 3.0), ylim=(-350,    0), ytick_step= 50),
    "case_6_IPyC_tangential_stress": dict(xlim=(0, 3.0), ylim=( -20,  180), ytick_step= 20),
    "case_6_SiC_tangential_stress":  dict(xlim=(0, 3.0), ylim=(-350,  100), ytick_step= 50),
    "case_7_IPyC_tangential_stress": dict(xlim=(0, 3.0), ylim=( -20,  200), ytick_step= 20),
    "case_7_SiC_tangential_stress":  dict(xlim=(0, 3.0), ylim=(-500,  100), ytick_step=100),

    # ── case 8 ──────────────────────────────────────────────────────────────
    "case_8_IPyC_tangential_stress": dict(xlim=(0, 3.0), ylim=( -50,  250), ytick_step= 50),
    "case_8_SiC_tangential_stress":  dict(xlim=(0, 3.0), ylim=(-700,  100), ytick_step=100),
}

# ---------------------------------------------------------------------------
# Per-code colour palette
#
# Canonical names are lower-case tokens; the helper get_series_style() maps
# any variant (different capitalisation, separators) to the right entry.
# Add a row whenever a new code appears in a benchmark figure.
# ---------------------------------------------------------------------------

import re as _re

# Canonical name → matplotlib kwargs
SERIES_STYLE = {
    "INL":    dict(color="#1f77b4", lw=1.5),   # US – Idaho National Laboratory
    "GA":     dict(color="#ff7f0e", lw=1.5),   # US – General Atomics
    "RF":     dict(color="#2ca02c", lw=1.5),   # Russian Federation code name
    "UK":     dict(color="#9467bd", lw=1.5),   # United Kingdom
    "Korea":  dict(color="#8c564b", lw=1.5),   # Korea
    "France": dict(color="#e377c2", lw=1.5),   # France
    "Japan":  dict(color="#7f7f7f", lw=1.5),   # Japan
    "Turkey": dict(color="#bcbd22", lw=1.5),   # Turkey
    "Russia": dict(color="#17becf", lw=1.5),   # Russia
}

# Aliases: known spelling/formatting variants → canonical name above.
# Case and separator variants are handled automatically by _norm(); only add
# entries here when the *normalised* strings still differ (e.g. "us inl" ≠ "inl").
_ALIASES = {
    "US INL":   "INL",
    "USA/INL":  "INL",
    "US-INL":   "INL",
    "U.S. INL": "INL",
    "US GA":    "GA",
    "USA/GA":   "GA",
    "US-GA":    "GA",
    "U.S. GA":  "GA",
    "U.K.":     "UK",   # normalises to "u k" ≠ "uk" without alias
}


def _norm(name: str) -> str:
    """Lower-case; collapse spaces, hyphens, slashes, dots to a single space."""
    return _re.sub(r"[\s\-/._]+", " ", name.strip()).lower()


_NORM_STYLE = {_norm(k): v for k, v in SERIES_STYLE.items()}
_NORM_ALIAS = {_norm(k): _norm(v) for k, v in _ALIASES.items()}


def get_series_style(name: str, fallback=None) -> dict:
    """Return matplotlib kwargs for *name*, resolving aliases and ignoring case/separators."""
    n = _norm(name)
    canonical = _NORM_ALIAS.get(n, n)
    return _NORM_STYLE.get(canonical, fallback if fallback is not None else dict(lw=1.5))


# OFFBEAT is always plotted as a thick red line on top of reference curves.
OFFBEAT_STYLE = dict(color="red", lw=2.5, zorder=5)
