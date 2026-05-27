"""
Shared pytest helpers for the IAEA HTGR FP diffusion benchmark.

Case-group mapping
------------------
  case_1*, case_2*  →  cases_1_2/
  case_3*           →  cases_3/
  case_4*           →  cases_4/
  case_5*           →  cases_5/

Tests invoke run.py via subprocess and then read postProcessing outputs.
"""

import os
import glob
import sys
import subprocess
import functools

import numpy as np
import pytest

_ROOT = os.path.dirname(os.path.abspath(__file__))

_GROUP = {
    "1": "cases_1_2",
    "2": "cases_1_2",
    "3": "cases_3",
    "4": "cases_4",
    "5": "cases_5",
}


def _cases_dir(case_id: str) -> str:
    """Return the group directory for a given case_id (e.g. 'case_3a' → 'cases_3')."""
    # case_id format: "case_Na..." where N is the benchmark group number
    num = case_id.split("_")[1][0]   # "case_3a" → "3"
    group = _GROUP.get(num)
    if group is None:
        raise ValueError(f"Cannot determine group directory for '{case_id}'")
    return os.path.join(_ROOT, group)


def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow (400 h FP diffusion simulations). "
        "Deselect with -m 'not slow'.",
    )


# ── Simulation runner ─────────────────────────────────────────────────────────

@functools.lru_cache(maxsize=None)
def run_case(case_id: str) -> str:
    """Run case_id via run.py in its group directory. Returns the case output path."""
    cwd = _cases_dir(case_id)
    subprocess.run([sys.executable, "run.py", case_id], cwd=cwd, check=True)
    return os.path.join(cwd, case_id)


# ── Result reader ─────────────────────────────────────────────────────────────

def _read_fp_dat(case_id: str) -> tuple[list[str], np.ndarray]:
    """Return (species_list, data_array) where data_array columns are [time, sp1, sp2, ...]."""
    case_dir = os.path.join(_cases_dir(case_id), case_id)
    pattern = os.path.join(case_dir, "postProcessing", "fpRelease", "*", "writeFpRelease.dat")
    files = glob.glob(pattern)
    if not files:
        pytest.skip(f"No fpRelease dat found for {case_id} — run the case first")
    files.sort(key=lambda p: float(os.path.basename(os.path.dirname(p))))

    species, rows = [], []
    with open(files[0]) as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                parts = s.lstrip("#").split()
                if parts and parts[0] == "Time":
                    species = parts[1:]
                continue
            parts = s.split()
            if len(parts) >= 2:
                rows.append([float(v) for v in parts])

    if not rows:
        pytest.skip(f"Empty fpRelease file for {case_id}")

    return species, np.array(rows)


def read_fp_release(case_id: str) -> dict[str, float]:
    """Return final fractional release per species {sp: value}."""
    species, arr = _read_fp_dat(case_id)
    return {sp: arr[-1, i + 1] for i, sp in enumerate(species)}


def read_fp_release_at(case_id: str, t: float) -> dict[str, float]:
    """Return fractional release per species at the time step closest to t (seconds)."""
    species, arr = _read_fp_dat(case_id)
    idx = int(np.argmin(np.abs(arr[:, 0] - t)))
    return {sp: arr[idx, i + 1] for i, sp in enumerate(species)}
