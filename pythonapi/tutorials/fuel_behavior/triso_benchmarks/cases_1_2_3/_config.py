"""
Geometry and material data for elastic TRISO benchmark cases 1–3.
No foamForNuclear dependency — imported by case.py and analytical.py.
"""

MAT = {
    "SiC": dict(
        density=3200.0, heatCapacity=1.0, conductivity=1.0, emissivity=0.0,
        YoungModulus=3.7e11, PoissonRatio=0.13, thermalExpansion=4.9e-6, Tref=293.0,
    ),
    "IPyC": dict(
        density=1900.0, heatCapacity=1.0, conductivity=1.0, emissivity=0.0,
        YoungModulus=3.96e10, PoissonRatio=0.33, thermalExpansion=5.5e-6, Tref=293.0,
    ),
}

CASES = {
    # case_1: SiC only — thermal gradient + pressure
    "case_1": {
        "inner_radius": 350.0,
        "layers": [("SiC", 385.0, 120)],
        "mat_names": ["SiC"],
        "T_int": 393.0, "T_inner": 793.0, "T_outer": 393.0,
        "p_inner": 25e6, "p_outer": 0.1e6,
        "endTime": 10.0, "deltaT": 1.0, "minDeltaT": 0.001,
    },
    # case_2: IPyC only — uniform temperature
    "case_2": {
        "inner_radius": 350.0,
        "layers": [("IPyC", 440.0, 80)],
        "mat_names": ["IPyC"],
        "T_int": 1273.0, "T_inner": 1273.0, "T_outer": 1273.0,
        "p_inner": 25e6, "p_outer": 0.1e6,
        "endTime": 10.0, "deltaT": 1.0, "minDeltaT": 0.001,
    },
    # case_3: IPyC + SiC — uniform temperature
    "case_3": {
        "inner_radius": 350.0,
        "layers": [("IPyC", 390.0, 80), ("SiC", 425.0, 80)],
        "mat_names": ["IPyC", "SiC"],
        "T_int": 1273.0, "T_inner": 1273.0, "T_outer": 1273.0,
        "p_inner": 25e6, "p_outer": 0.1e6,
        "endTime": 10.0, "deltaT": 1.0, "minDeltaT": 0.001,
    },
}
