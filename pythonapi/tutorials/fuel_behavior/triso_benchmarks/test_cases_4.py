"""
Verification tests for TRISO PyC benchmark cases 4a–4d.

Asserts that the final-time tangential stress at the inner face of IPyC and
SiC falls within the CRP-6 code-comparison range (±0.5 MPa tolerance).

Cases run for 1000 days (irradiation) and are marked @pytest.mark.slow.
Deselect with:  pytest -m "not slow"
"""

import pytest
from conftest import read_last_sigma, sigma_at, run_case

um  = 1e-6
TOL = 0.5   # MPa — accounts for floating-point proximity to range boundaries

# IAEA TECDOC-1674 §9 CRP-6 code-comparison ranges [MPa]
CRP6 = {
    ("case_4a", "IPyC"): ( 925,  970),
    ("case_4a", "SiC"):  (-850, -775),
    ("case_4b", "IPyC"): ( -25,  -25),
    ("case_4b", "SiC"):  ( 138,  142),
    ("case_4c", "IPyC"): (  25,   27),
    ("case_4c", "SiC"):  (  83,   92),
    ("case_4d", "IPyC"): (  25,   35),
    ("case_4d", "SiC"):  (  71,   88),
}

# Inner face radii: inner_radius=350 µm, IPyC thickness=40 µm
LAYER_R = {
    "IPyC": 350 * um,
    "SiC":  390 * um,
}


@pytest.mark.slow
@pytest.mark.parametrize("case_id,layer", [
    ("case_4a", "IPyC"), ("case_4a", "SiC"),
    ("case_4b", "IPyC"), ("case_4b", "SiC"),
    ("case_4c", "IPyC"), ("case_4c", "SiC"),
    ("case_4d", "IPyC"), ("case_4d", "SiC"),
])
def test_tangential_stress(case_id, layer):
    case_dir = run_case("cases_4", case_id)
    r, _, st = read_last_sigma(case_dir)
    val = sigma_at(r, st, LAYER_R[layer])
    lo, hi = CRP6[(case_id, layer)]
    assert lo - TOL <= val <= hi + TOL, (
        f"{case_id} {layer} σ_θ = {val:.1f} MPa, CRP-6 range = [{lo}, {hi}] MPa"
    )
