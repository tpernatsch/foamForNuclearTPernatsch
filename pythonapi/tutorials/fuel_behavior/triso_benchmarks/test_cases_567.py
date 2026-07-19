"""
Verification tests for TRISO PyC benchmark cases 5–7.

Asserts that the final-time tangential stress at the inner face of IPyC and
SiC falls within the CRP-6 code-comparison range (±0.5 MPa tolerance).

Cases run for 1000 days (irradiation) and are marked @pytest.mark.slow.
Deselect with:  pytest -m "not slow"
"""

import pytest
from conftest import read_last_sigma, sigma_at, run_case

um  = 1e-6
TOL = 0.5   # MPa

# IAEA TECDOC-1674 §9 CRP-6 code-comparison ranges [MPa]
CRP6 = {
    ("case_5", "IPyC"): ( 40,  58),
    ("case_5", "SiC"):  (-56, -28),
    ("case_6", "IPyC"): ( 27,  38),
    ("case_6", "SiC"):  ( 28,  48),
    ("case_7", "IPyC"): ( 37,  50),
    ("case_7", "SiC"):  ( 10,  25),
}

# Inner-face radii vary by case (case_5 has smaller inner_radius=275 µm).
# IPyC thickness=40 µm, SiC thickness=35 µm for all cases.
LAYER_R = {
    "case_5": {"IPyC": 275 * um, "SiC": 315 * um},
    "case_6": {"IPyC": 350 * um, "SiC": 390 * um},
    "case_7": {"IPyC": 350 * um, "SiC": 390 * um},
}


@pytest.mark.slow
@pytest.mark.parametrize("case_id,layer", [
    ("case_5", "IPyC"), ("case_5", "SiC"),
    ("case_6", "IPyC"), ("case_6", "SiC"),
    ("case_7", "IPyC"), ("case_7", "SiC"),
])
def test_tangential_stress(case_id, layer):
    case_dir = run_case("cases_5_6_7", case_id)
    r, _, st = read_last_sigma(case_dir)
    val = sigma_at(r, st, LAYER_R[case_id][layer])
    lo, hi = CRP6[(case_id, layer)]
    assert lo - TOL <= val <= hi + TOL, (
        f"{case_id} {layer} σ_θ = {val:.1f} MPa, CRP-6 range = [{lo}, {hi}] MPa"
    )
