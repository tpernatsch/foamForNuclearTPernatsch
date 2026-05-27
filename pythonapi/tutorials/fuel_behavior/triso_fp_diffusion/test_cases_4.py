"""
Verification tests for IAEA HTGR FP diffusion benchmark cases 4a–4d.

Cases 4 mirror cases 3 but add pre-irradiation fast fluence
Γ = 2×10²¹ n/cm² (Table 10.5), which enhances D_Cs(SiC) via exp(Γ × 2.2×10⁻²²).

Cases 4a–4c track Cs-137 and Ag-110m.
Case 4d: SiC cracks at 1600°C (t=0), IPyC+OPyC crack at 1800°C (t_switch)
→ tracks Cs-137 and Kr-85 (same failure sequence and species as case 3e).

Reference values from Table 10.9.
"""

import pytest
from conftest import run_case, read_fp_release

# fmt: off
@pytest.mark.parametrize("case_id,species,ref,tol_rel,tol_abs", [
    # Wide Cs tolerance: Table 10.9 shows 4a spread of 1.64e-4 – 1.47e-3
    ("case_4a", "Cs", 4.10e-4, 1.00, None),
    ("case_4a", "Ag", 0.43,    0.30, None),
    ("case_4b", "Cs", 0.21,    0.10, None),
    ("case_4b", "Ag", 0.89,    0.10, None),
    pytest.param("case_4c", "Cs", 0.23, 0.10, None,  marks=pytest.mark.slow),
    pytest.param("case_4c", "Ag", 0.93, 0.10, None,  marks=pytest.mark.slow),
    # case_4d: Kr-85 (not Ag), same failure sequence as case_3e
    pytest.param("case_4d", "Cs", 1.0,  None, 0.05,  marks=pytest.mark.slow),
    pytest.param("case_4d", "Kr", 1.0,  None, 0.05,  marks=pytest.mark.slow),
])
# fmt: on
def test_fp_release(case_id, species, ref, tol_rel, tol_abs):
    run_case(case_id)
    final = read_fp_release(case_id)[species]

    if tol_abs is not None:
        assert abs(final - ref) <= tol_abs, (
            f"{case_id} {species} = {final:.6g}, expected ≈ {ref:.4g}  "
            f"(abs err = {abs(final - ref):.4g}, tol = {tol_abs})"
        )
    else:
        rel_err = abs(final - ref) / ref
        assert rel_err < tol_rel, (
            f"{case_id} {species} = {final:.6g}, expected {ref:.4g}  "
            f"(rel err = {rel_err * 100:+.1f} %, tol = {tol_rel * 100:.0f} %)"
        )
