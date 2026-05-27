"""
IAEA HTGR FP diffusion benchmark — cases 4a–4d (pre-irradiated TRISO).

Irradiation conditions (Table 10.5): 500 efpd at 1000°C, 10% FIMA,
fast fluence Γ = 2×10²⁵ n/m² = 2×10²¹ n/cm².

Each simulation covers the full history in two phases:
  Phase 1 — irradiation at 1000°C for 500 efpd, Cs and Ag produced from fission.
  Phase 2 — heating test (step change to target temperature, no fission power).

Cases 4a–4c track Cs-137 and Ag-110m.
Case 4d (same failure sequence as case 3e): Cs-137 and Kr-85.

Reference values from Table 10.9.
"""

# ── Kernel material parameters (must match materials.UO2 in case.py) ──────
_RHO_TD         = 10960.0   # kg/m³  theoretical density
_DENSITY_FRAC   = 0.95      # density fraction
_M_UO2          = 0.270     # kg/mol
_N_A            = 6.022e23  # /mol
_E_FISSION      = 3.2e-11   # J per fission

# ── Irradiation phase ─────────────────────────────────────────────────────
T_IRRAD         = 1273.15          # 1000°C
EFPD_IRRAD      = 500
T_IRRAD_END     = EFPD_IRRAD * 86400  # 43 200 000 s

# Number density of U atoms (= UO2 molecule density)
_N_f = _RHO_TD * _DENSITY_FRAC * _N_A / _M_UO2   # /m³

# VHGR [W/m³] for 10% FIMA in T_IRRAD_END seconds
VHGR = 0.10 * _N_f * _E_FISSION / T_IRRAD_END

# Fast flux rate [n/cm²/s] to accumulate Γ = 2×10²¹ n/cm² in T_IRRAD_END s
FAST_FLUX_RATE = 2.0e21 / T_IRRAD_END   # n/cm²/s

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
    # ── Single-phase heating ───────────────────────────────────────────────
    "case_4a": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T_heat":    (1873.15,),           # step to 1600°C at end of irradiation
        "endTime":   T_IRRAD_END + 720000, # 500 efpd + 200 h
        "deltaT":    3600,
        # Table 10.9 — wide participant spread (Cs: 1.64e-4 – 1.47e-3, Ag: 0.27 – 0.55)
        "expected":  {"Cs": 4.10e-4, "Ag": 0.43},   # USA/INL representative
    },
    "case_4b": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T_heat":    (2073.15,),           # step to 1800°C
        "endTime":   T_IRRAD_END + 720000,
        "deltaT":    3600,
        # Table 10.9 — tight Cs cluster; Ag: France outlier (0.58), others 0.87–0.95
        "expected":  {"Cs": 0.21, "Ag": 0.89},       # cluster centres
    },
    # ── Two-phase heating ─────────────────────────────────────────────────
    "case_4c": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T_heat":    (1873.15, 2073.15),   # step to 1600°C, then step to 1800°C
        "endTime":   T_IRRAD_END + 1440000,
        "deltaT":    3600,
        # Table 10.9 — tight clusters
        "expected":  {"Cs": 0.23, "Ag": 0.93},
    },
    "case_4d": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "T_heat":    (1873.15, 2073.15),
        "endTime":   T_IRRAD_END + 1440000,
        "deltaT":    3600,
        # Table 10.5: SiC cracks at 1600°C (= start of heating = T_IRRAD_END),
        # IPyC+OPyC crack at 1800°C (= temperature step = T_IRRAD_END + 200 h)
        "sic_failure_time": T_IRRAD_END,
        "pyc_failure_time": T_IRRAD_END + 720000,
        "include_kr": True,              # Table 10.5: Kr-85 tracked (same as case_3e)
        # Table 10.9: all coatings fail → Cs → 1, Kr → 1
        "expected":  {"Cs": 1.0, "Kr": 1.0},
    },
}
