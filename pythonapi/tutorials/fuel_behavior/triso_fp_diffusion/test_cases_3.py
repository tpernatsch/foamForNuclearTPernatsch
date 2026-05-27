"""
Verification tests for IAEA HTGR FP diffusion benchmark cases 3a–3e.

Cases 3a/3b (Table 10.8) — full TRISO, single-phase heating, 200 h:
  3a (1600°C): participant spread 6.59×10⁻⁵ – 1.15×10⁻³ → 50 % relative.
  3b (1800°C): cluster 0.203–0.218 → 10 % relative.

Cases 3c–3e — two-phase heating, 400 h total:
  3c: cluster 0.220–0.239 → 10 % relative.
  3d/3e: Cs → 1 once coatings fail → ±0.05 absolute.
"""

import pytest
from conftest import run_case, read_fp_release

EXPECTED_CS = {
    "case_3a": 1.54447e-4,
    "case_3b": 0.210283,
    "case_3c": 0.225066,
    "case_3d": 1.0,
    "case_3e": 1.0,
}


@pytest.mark.parametrize("case_id,tol_rel,tol_abs", [
    ("case_3a", 0.50,  None  ),
    ("case_3b", 0.10,  None  ),
    pytest.param("case_3c", 0.10, None,  marks=pytest.mark.slow),
    pytest.param("case_3d", None, 0.05,  marks=pytest.mark.slow),
    pytest.param("case_3e", None, 0.05,  marks=pytest.mark.slow),
])
def test_cs_release(case_id, tol_rel, tol_abs):
    run_case(case_id)
    final = read_fp_release(case_id)["Cs"]
    ref = EXPECTED_CS[case_id]

    if tol_abs is not None:
        assert abs(final - ref) <= tol_abs, (
            f"{case_id} Cs = {final:.6g}, expected ≈ {ref:.4g}  "
            f"(abs err = {abs(final-ref):.4g}, tol = {tol_abs})"
        )
    else:
        rel_err = abs(final - ref) / ref
        assert rel_err < tol_rel, (
            f"{case_id} Cs = {final:.6g}, expected {ref:.4g}  "
            f"(rel err = {rel_err*100:+.1f} %, tol = {tol_rel*100:.0f} %)"
        )
