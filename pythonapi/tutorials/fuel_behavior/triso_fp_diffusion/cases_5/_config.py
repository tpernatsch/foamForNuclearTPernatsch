"""
IAEA HTGR FP diffusion benchmark — case 5 (full TRISO, Cs + Ag, cyclic irradiation).

Irradiation conditions (Table 10.5):
  10 cycles × 100 efpd, temperature ramp 600→1000°C per cycle,
  10 % FIMA, fast fluence Γ = 2×10²⁵ n/m² = 2×10²¹ n/cm².

Case 5a: results reported immediately after irradiation (t = 24 000 h).
Case 5b: irradiation + 200 h post-irradiation heating at 1600°C.

Reference values from Table 10.10.
"""

# ── Kernel material parameters (must match materials.UO2 in case.py) ──────
_RHO_TD         = 10960.0   # kg/m³  theoretical density
_DENSITY_FRAC   = 0.95      # density fraction
_M_UO2          = 0.270     # kg/mol
_N_A            = 6.022e23  # /mol
_E_FISSION      = 3.2e-11   # J per fission

# ── Cyclic irradiation phase ──────────────────────────────────────────────
N_CYCLES        = 10
EFPD_PER_CYCLE  = 100
CYCLE_DURATION  = EFPD_PER_CYCLE * 86400    # 8 640 000 s  (100 efpd)
T_IRRAD_END     = N_CYCLES * CYCLE_DURATION  # 86 400 000 s (24 000 h)

T_IRRAD_LOW     = 873.15     # 600°C  — start of each irradiation cycle
T_IRRAD_HIGH    = 1273.15    # 1000°C — end of each irradiation cycle
T_HEATING       = 1873.15    # 1600°C — post-irradiation heating (case_5b only)

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
    # ── Single simulation covering full history: irradiation + heating ─────────
    # Case 5a (irradiation-end snapshot) and case 5b (post-heating) are both
    # extracted from this one run in post-processing.
    # deltaT = 1 h = 3 600 s  →  24 000 steps irradiation + 200 steps heating.
    "case_5b": {
        "layers":    _FULL_TRISO,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC", "SiC"],
        "t_heating": T_IRRAD_END,                # heating phase starts here
        "endTime":   T_IRRAD_END + 200 * 3600,   # 87 120 000 s (24 200 h)
        "deltaT":    3600,                        # 1 h per step throughout
        # Table 10.10:
        #   5a (at T_IRRAD_END): Ag ~2×10⁻⁵ (Cs reference unreliable)
        #   5b (at end of heating): Cs ~6.5×10⁻⁴, Ag ~0.42
        "expected_5a": {"Ag": 2e-5},
        "expected":    {"Cs": 6.5e-4, "Ag": 0.42},
    },
}
