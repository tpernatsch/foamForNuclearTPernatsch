"""
Verification tests for IAEA HTGR FP diffusion benchmark case 5 (Cs + Ag).

A single simulation (case_5b) covers the full history:
  Phase 1 — 10 cycles × 100 efpd, T ramp 600→1000°C, Γ = 2×10²¹ n/cm²
  Phase 2 — 200 h heating at 1600°C (step change)

The case_5a snapshot (end-of-irradiation release) is read from case_5b at
T_IRRAD_END rather than from a separate simulation.

Reference values from Table 10.10.
Note: case_5a Cs has a 7-decade participant spread — no reliable reference;
that assertion is omitted.
"""

import pytest
from conftest import run_case, read_fp_release, read_fp_release_at
from cases_5._config import T_IRRAD_END


# ── case_5a snapshot (end of irradiation, read from case_5b output) ──────────
@pytest.mark.slow
def test_5a_ag_end_of_irradiation():
    run_case("case_5b")
    val = read_fp_release_at("case_5b", T_IRRAD_END)["Ag"]
    ref, tol_rel = 2e-5, 1.00
    rel_err = abs(val - ref) / ref
    assert rel_err < tol_rel, (
        f"case_5a (snapshot) Ag = {val:.6g}, expected ≈ {ref:.4g}  "
        f"(rel err = {rel_err * 100:+.1f} %, tol = {tol_rel * 100:.0f} %)"
    )


# ── case_5b final values (end of heating) ────────────────────────────────────
# fmt: off
@pytest.mark.parametrize("species,ref,tol_rel", [
    ("Cs", 6.5e-4, 0.80),
    ("Ag", 0.42,   0.50),
])
# fmt: on
@pytest.mark.slow
def test_5b_fp_release(species, ref, tol_rel):
    run_case("case_5b")
    final = read_fp_release("case_5b")[species]
    rel_err = abs(final - ref) / ref
    assert rel_err < tol_rel, (
        f"case_5b {species} = {final:.6g}, expected {ref:.4g}  "
        f"(rel err = {rel_err * 100:+.1f} %, tol = {tol_rel * 100:.0f} %)"
    )
