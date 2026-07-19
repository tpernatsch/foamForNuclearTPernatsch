"""
Geometry, material and irradiation data for PyC TRISO benchmark cases 4a–4d.
No foamForNuclear dependency — imported by case.py.

All cases: IPyC/SiC two-layer geometry, T=1273 K uniform, 1000-day irradiation.
  4a: swelling only (constant isotropic rate), no creep
  4b: creep only (constant), no swelling
  4c: constant swelling + constant creep
  4d: dose-dependent anisotropic swelling (correlation a) + constant creep
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
}

# Swelling correlation (a) — same as cases 5 and 6
PYC_SWELLING_COEFFICIENTS = {
    "A": {
        "radial":     [-2.22642e-2,  2.00861e-2, -7.77024e-3,  1.36334e-3, 0.0, 0.0],
        "tangential": [-1.91253e-2,  2.63307e-3,  1.69251e-3, -3.53804e-4, 0.0, 0.0],
    },
}

# Fast flux: constant 3.4722e13 n/cm²/s over 1000 days
PYC_FAST_FLUX = dict(
    timePoints=[0.0, 8.64e7],
    fastFlux=[3.4722e13, 3.4722e13],
    materials=["IPyC", "SiC"],
)

# Layer thicknesses in μm: IPyC=40, SiC=35 (no OPyC)
_LAYERS    = [("IPyC", 40.0, 60), ("SiC", 35.0, 60)]
_MAT_NAMES = ["IPyC", "SiC"]
_TIME      = dict(endTime=8.64e7, deltaT=8.64e3, minDeltaT=1.0, maxDeltaT=4.32e5)

CASES = {
    # case_4a: swelling only, constant isotropic rate -0.005 per 10²⁵ n/m²
    "case_4a": dict(
        inner_radius=350.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        T_int=1273.0, T_inner=1273.0, T_outer=1273.0,
        p_inner=25e6, p_outer=0.1e6,
        pyc_Tref=293.0,
        pyc_swelling="constant", pyc_swelling_rate=-0.005,
        pyc_creep=None,
        **_TIME,
    ),
    # case_4b: creep only, constant coefficient 2.71e-4 MPa/(10²⁵ n/m²)
    "case_4b": dict(
        inner_radius=350.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        T_int=1273.0, T_inner=1273.0, T_outer=1273.0,
        p_inner=25e6, p_outer=0.1e6,
        pyc_Tref=293.0,
        pyc_swelling=None,
        pyc_creep="constant", pyc_creep_coefficient=2.71e-4,
        **_TIME,
    ),
    # case_4c: constant swelling + constant creep
    "case_4c": dict(
        inner_radius=350.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        T_int=1273.0, T_inner=1273.0, T_outer=1273.0,
        p_inner=25e6, p_outer=0.1e6,
        pyc_Tref=293.0,
        pyc_swelling="constant", pyc_swelling_rate=-0.005,
        pyc_creep="constant", pyc_creep_coefficient=2.71e-4,
        **_TIME,
    ),
    # case_4d: anisotropic dose-dependent swelling (correlation a) + constant creep
    "case_4d": dict(
        inner_radius=350.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        T_int=1273.0, T_inner=1273.0, T_outer=1273.0,
        p_inner=25e6, p_outer=0.1e6,
        pyc_Tref=293.0,
        pyc_swelling="correlation", pyc_swelling_set="A",
        pyc_creep="constant", pyc_creep_coefficient=2.71e-4,
        **_TIME,
    ),
}
