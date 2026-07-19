"""
Smoke test for TRISO PyC benchmark case 8 (10-cycle temperature profile).

No CRP-6 numerical range is available yet — the benchmark comparison in
TECDOC-1674 is graphical only.  This test asserts the simulation ran and
produced output.  Promote to a range assertion once values are digitized
from the reference figure with WebPlotDigitizer.

Marked @pytest.mark.slow (1000-day irradiation with daily time steps).
Deselect with:  pytest -m "not slow"
"""

import os
import glob
import pytest
from conftest import read_last_sigma, run_case


@pytest.mark.slow
def test_output_exists():
    case_dir = run_case("cases_8", "case_8")
    pattern = os.path.join(
        case_dir, "postProcessing", "radialProfile", "*", "line_T_sigma.csv"
    )
    files = glob.glob(pattern)
    assert files, f"No radialProfile CSV produced in {case_dir}"


@pytest.mark.slow
def test_output_non_empty():
    case_dir = run_case("cases_8", "case_8")
    r, sr, st = read_last_sigma(case_dir)
    assert len(r) > 0, "radialProfile CSV is empty"
    assert len(sr) == len(r)
    assert len(st) == len(r)
