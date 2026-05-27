"""
IAEA HTGR FP diffusion benchmark — cases 1a/1b (bare kernel) and 2a/2b (kernel+Buffer+IPyC).
No pre-irradiation fluence; fixed temperature throughout.

Case 1  (Table 10.6) — bare kernel, 200 h:
  Analytical solution:  1a = 0.4673  (1200°C),  1b = 0.99999959  (1600°C).

Case 2  (Table 10.7) — kernel + Buffer + IPyC, 200 h:
  Participant cluster:  2a 0.026–0.030  (1200°C),  2b 0.968–0.996  (1600°C).
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
}

_KERNEL_ONLY = [("Kernel", 250.0, 50)]

_KERNEL_BUFFER_IPYC = [
    ("Kernel", 250.0, 50),
    ("Buffer", 350.0, 30),
    ("IPyC",   390.0, 15),
]

CASES = {
    # ── Cases 1 : bare kernel ─────────────────────────────────────────────────
    "case_1a": {
        "layers":    _KERNEL_ONLY,
        "mat_names": ["Kernel"],
        "T":         1473.15,    # 1200°C
        "fluence":   0.0,
        "endTime":   720000,     # 200 h
        "deltaT":    3600,
        "expected":  {"Cs": 0.4673},       # analytical solution (Table 10.6)
    },
    "case_1b": {
        "layers":    _KERNEL_ONLY,
        "mat_names": ["Kernel"],
        "T":         1873.15,    # 1600°C
        "fluence":   0.0,
        "endTime":   720000,
        "deltaT":    3600,
        "expected":  {"Cs": 0.99999959},   # analytical solution (Table 10.6)
    },
    # ── Cases 2 : kernel + Buffer + IPyC ─────────────────────────────────────
    "case_2a": {
        "layers":    _KERNEL_BUFFER_IPYC,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC"],
        "T":         1473.15,    # 1200°C
        "fluence":   0.0,
        "endTime":   720000,
        "deltaT":    3600,
        "expected":  {"Cs": 0.026},        # Table 10.7 cluster 0.026–0.030
    },
    "case_2b": {
        "layers":    _KERNEL_BUFFER_IPYC,
        "mat_names": ["Kernel", "Buffer", "IPyC|OPyC"],
        "T":         1873.15,    # 1600°C
        "fluence":   0.0,
        "endTime":   720000,
        "deltaT":    3600,
        "expected":  {"Cs": 0.994714},     # Table 10.7 cluster 0.968–0.996
    },
}
