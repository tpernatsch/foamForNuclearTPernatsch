"""
Verification tests for IAEA HTGR FP diffusion benchmark cases 1 & 2.

Case 1  (Table 10.6) — bare kernel, 200 h:
  Analytical solution: 1a = 0.4673, 1b = 0.99999959.
  → 5 % relative (1a) / ±0.03 absolute (1b).

Case 2  (Table 10.7) — kernel + Buffer + IPyC, 200 h:
  Participant cluster: 2a 0.026–0.030, 2b 0.968–0.996.
  (Harmonic Laplacian scheme required for 2a.)
  → 10 % relative.
"""

import pytest
from conftest import run_case, read_fp_release

EXPECTED_CS = {
    "case_1a": 0.4673,
    "case_1b": 0.99999959,
    "case_2a": 0.026,
    "case_2b": 0.994714,
}


@pytest.mark.parametrize("case_id,tol_rel,tol_abs", [
    ("case_1a", 0.05,  None  ),
    ("case_1b", None,  0.03  ),
    ("case_2a", 0.10,  None  ),
    ("case_2b", 0.05,  None  ),
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
