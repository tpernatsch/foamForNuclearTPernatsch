"""
Geometry, material and cycling data for PyC TRISO benchmark case 8.
No foamForNuclear dependency — imported by case.py.

Case 8: 10 temperature cycles over 1000 days.
  - Temperature ramps linearly from 873 K to 1273 K over each 100-day cycle,
    then drops immediately back to 873 K (IAEA TECDOC-1674 §9).
  - Inner pressure and fast fluence follow a sawtooth profile given in the TECDOC.
"""

MAT = {
    "IPyC": dict(
        density=1900.0, heatCapacity=1.0, conductivity=1.0, emissivity=0.0,
        YoungModulus=3.96e10, PoissonRatio=0.33, thermalExpansion=5.5e-6, Tref=293.0,
    ),
    "SiC": dict(
        density=3200.0, heatCapacity=1.0, conductivity=1.0, emissivity=0.0,
        YoungModulus=3.7e11, PoissonRatio=0.13, thermalExpansion=4.9e-6, Tref=293.0,
    ),
    "OPyC": dict(
        density=1900.0, heatCapacity=1.0, conductivity=1.0, emissivity=0.0,
        YoungModulus=3.96e10, PoissonRatio=0.33, thermalExpansion=5.5e-6, Tref=293.0,
    ),
}

PYC_SWELLING_COEFFICIENTS = {
    "C": {
        "radial":     [-1.80613e-2,  9.82884e-3, -2.25937e-3,  4.03266e-4, 0.0, 0.0],
        "tangential": [-1.78392e-2,  1.71315e-3,  2.32979e-3, -4.91648e-4, 0.0, 0.0],
    },
}

# Fast flux: constant 3.4722e13 n/cm²/s over 1000 days
PYC_FAST_FLUX = dict(
    timePoints=[0.0, 8.64e7],
    fastFlux=[3.4722e13, 3.4722e13],
    materials=["IPyC", "SiC", "OPyC"],
)

# Layer thicknesses in μm: IPyC=40, SiC=35, OPyC=40
_LAYERS    = [("IPyC", 40.0, 8), ("SiC", 35.0, 8), ("OPyC", 40.0, 8)]
_MAT_NAMES = ["IPyC", "SiC", "OPyC"]

# ------------------------------------------------------------------ #
# Cycling profile (fluence in 10²⁵ n/m², converted to time in s)
# ------------------------------------------------------------------ #
_ff = PYC_FAST_FLUX["fastFlux"][0]   # n/cm²/s

_PHI = [
    0,    0.29,
    0.30, 0.59,
    0.60, 0.89,
    0.90, 1.19,
    1.20, 1.49,
    1.50, 1.79,
    1.80, 2.09,
    2.10, 2.39,
    2.40, 2.69,
    2.70, 2.99,
    3.00,
]

_TIMES = [phi * 1e25 / (_ff * 1e4) for phi in _PHI]   # seconds

# Temperature: 873 K at cycle start (even index), 1273 K at ramp end (odd index)
_T_PROFILE = [873.0 if i % 2 == 0 else 1273.0 for i in range(len(_PHI))]

# Inner pressure [MPa] at each profile point, converted to Pa
_P_MPa = [
     0.,    0.14,
     0.02,  0.94,
     0.04,  2.59,
     0.07,  4.87,
     0.10,  7.64,
     0.14, 10.79,
     0.20, 14.26,
     0.26, 17.99,
     0.33, 21.96,
     0.41, 26.13,
     0.50,
]
_P_PROFILE = [p * 1e6 for p in _P_MPa]   # Pa

CASES = {
    "case_8": dict(
        inner_radius=350.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        T_list=list(zip(_TIMES, _T_PROFILE)),
        p_inner_list=list(zip(_TIMES, _P_PROFILE)),
        p_outer=0.1e6,
        pyc_Tref=298.15, pyc_swelling_set="C",
        pyc_creep_coefficients=[4.386e-4, -9.70e-7, 8.0294e-10],
        endTime=8.64e7, deltaT=8.64e3, minDeltaT=1.0, maxDeltaT=8.64e4,
    ),
}
