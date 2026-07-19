"""
Verification tests for TRISO elastic benchmark cases 1–3.

Compares σ_θ at the inner face of each layer against the TECDOC-1674 §9.4
Lamé analytical solution (2 % relative tolerance).  Case 3 also checks
σ_r at the IPyC/SiC interface.

These cases are purely elastic (no swelling/creep), so they run in seconds
and carry no @pytest.mark.slow marker.
"""

import pytest
from conftest import read_last_sigma, sigma_at, run_case

um      = 1e-6
TOL_REL = 0.02   # 2 % relative error

# TECDOC-1674 §9.4 reference tangential stresses [MPa] at r_target [µm]
REFS = {
    ("case_1", "SiC",  350): 125.19,
    ("case_2", "IPyC", 350):  50.20,
    ("case_3", "IPyC", 350):   8.80,
    ("case_3", "SiC",  390): 104.40,
}

# In cases_1_2_3 the layer tuple stores outer radii [µm], so IPyC/SiC
# interface for case_3 is at the outer radius of IPyC = 390 µm.
IFACE_R_CASE3      = 390 * um
IFACE_RADIAL_REF   = -18.8   # MPa  (σ_r = −Π, compressive contact pressure)


@pytest.mark.parametrize("case_id,layer,r_um", [
    ("case_1", "SiC",  350),
    ("case_2", "IPyC", 350),
    ("case_3", "IPyC", 350),
    ("case_3", "SiC",  390),
])
def test_tangential_stress(case_id, layer, r_um):
    case_dir = run_case("cases_1_2_3", case_id)
    r, _, st = read_last_sigma(case_dir)
    val = sigma_at(r, st, r_um * um)
    ref = REFS[(case_id, layer, r_um)]
    assert abs(val - ref) / abs(ref) < TOL_REL, (
        f"{case_id} {layer} σ_θ = {val:.2f} MPa, TECDOC ref = {ref:.2f} MPa"
    )


def test_interface_radial_case3():
    case_dir = run_case("cases_1_2_3", "case_3")
    r, sr, _ = read_last_sigma(case_dir)
    val = sigma_at(r, sr, IFACE_R_CASE3)
    ref = IFACE_RADIAL_REF
    assert abs(val - ref) / abs(ref) < TOL_REL, (
        f"case_3 IPyC/SiC interface σ_r = {val:.2f} MPa, TECDOC ref = {ref:.2f} MPa"
    )
