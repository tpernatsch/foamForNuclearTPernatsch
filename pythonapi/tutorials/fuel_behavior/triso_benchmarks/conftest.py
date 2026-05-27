"""
Shared helpers and fixtures for the TRISO benchmark verification suite.

Tests do NOT import from the case-group modules — they invoke run.py via
subprocess and then read the CSV outputs.  This avoids sys.path conflicts
between case groups that each have their own _config.py / case.py.
"""

import os
import glob
import sys
import subprocess
import functools

import numpy as np
import pandas as pd
import pytest

BENCH_DIR = os.path.dirname(os.path.abspath(__file__))
MPa = 1e6


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow (1000-day irradiation simulations). "
        "Deselect with -m 'not slow'.",
    )


# ── I/O helpers ─────────────────────────────────────────────────────────────

def read_last_sigma(case_dir):
    """Return (r, sigma_r, sigma_theta) numpy arrays [SI] from the last CSV."""
    pattern = os.path.join(
        case_dir, "postProcessing", "radialProfile", "*", "line_T_sigma.csv"
    )
    files = sorted(
        glob.glob(pattern),
        key=lambda p: float(os.path.basename(os.path.dirname(p)))
    )
    if not files:
        pytest.skip(f"No radialProfile CSV in {case_dir} — run the case first")
    df = pd.read_csv(files[-1])
    return df["x"].to_numpy(), df["sigma_0"].to_numpy(), df["sigma_3"].to_numpy()


def sigma_at(r_sim, sigma_arr, r_target):
    """Stress value [MPa] at the point nearest to r_target [m]."""
    return sigma_arr[np.argmin(np.abs(r_sim - r_target))] / MPa


# ── Simulation runner (cached so each case runs at most once per session) ───

@functools.lru_cache(maxsize=None)
def run_case(case_group, case_id):
    """
    Run a benchmark case via its run.py and return the case output directory.
    Results are cached — calling this twice for the same case is free.
    """
    group_dir = os.path.join(BENCH_DIR, case_group)
    subprocess.run(
        [sys.executable, "run.py", case_id],
        cwd=group_dir,
        check=True,
    )
    return os.path.join(group_dir, case_id)
