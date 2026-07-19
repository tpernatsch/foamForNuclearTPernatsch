"""
IAEA HTGR FP diffusion benchmark — cases 3a–3e (full TRISO, no pre-irradiation).

Cases 3a/3b — single-phase heating (200 h at 1600°C / 1800°C).
Cases 3c/3d/3e — two-phase heating: 200 h at 1600°C then 200 h at 1800°C (400 h total).
  3d: SiC cracks at the 200 h temperature step.
  3e: SiC cracked from t=0; IPyC+OPyC also crack at 200 h.

Reference: Table 10.8 of the IAEA HTGR FP diffusion benchmark.
"""

MAT = {
    "Buffer": dict(
        density=1880, conductivity=4.0, heatCapacity=720, emissivity=1.0,
        YoungModulus=3.96e10, PoissonRatio=0.33, thermalExpansion=5.5e-6, Tref=1608,
    ),
    "IPyC|OPyC": dict(
        density=1880, conductivity=4.0, heatCapacity=720, emissivity=1.0,
        YoungModulus=3.96e10, PoissonRatio=0.33, thermalExpansion=5.5e-6, Tref=1608,
    ),
    "SiC": dict(
        density=3200, conductivity=13.9, heatCapacity=620, emissivity=0.0,
        YoungModulus=3.7e11, PoissonRatio=0.13, thermalExpansion=4.9e-6, Tref=1608,
    ),
}

_FULL_TRISO = [
    ("Kernel", 250.0, 50),
    ("Buffer", 350.0, 30),
    ("IPyC",   390.0, 15),
    ("SiC",    425.0, 15),
    ("OPyC",   465.0, 15),
]

CASES = {
    # ── Single-phase heating ───────────────────────────────────────────────────
    "case_3a": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T":         1873.15,    # 1600°C
        "fluence":   0.0,
        "endTime":   720000,     # 200 h
        "deltaT":    3600,
        "expected":  {"Cs": 1.54447e-4},   # Table 10.8
    },
    "case_3b": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T":         2073.15,    # 1800°C
        "fluence":   0.0,
        "endTime":   720000,
        "deltaT":    3600,
        "expected":  {"Cs": 0.210283},     # Table 10.8 cluster 0.203–0.218
    },
    # ── Two-phase heating ─────────────────────────────────────────────────────
    "case_3c": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T":         (1873.15, 2073.15),   # 1600°C → 1800°C
        "t_switch":  720000,               # step at 200 h
        "fluence":   0.0,
        "endTime":   1440000,              # 400 h total
        "deltaT":    3600,
        "expected":  {"Cs": 0.225066},    # Table 10.8 cluster 0.220–0.239
    },
    "case_3d": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T":         (1873.15, 2073.15),
        "t_switch":  720000,
        "fluence":   0.0,
        "endTime":   1440000,
        "deltaT":    3600,
        "sic_failure_time": 720000,        # SiC cracks at 200 h
        "expected":  {"Cs": 1.0},          # Table 10.8: 0.999–1.000
    },
    "case_3e": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T":         (1873.15, 2073.15),
        "t_switch":  720000,
        "fluence":   0.0,
        "endTime":   1440000,
        "deltaT":    3600,
        "sic_failure_time": 0,             # SiC cracked from t=0 (when T=1600°C)
        "pyc_failure_time": 720000,        # IPyC+OPyC crack at 200 h (when T=1800°C)
        "include_kr": True,                # Table 10.5: tracks Cs-137 + Kr-85
        "expected":  {"Cs": 1.0, "Kr": 1.0},  # all coatings fail → both → 1
    },
}
