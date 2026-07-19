"""
Geometry, material and swelling data for PyC TRISO benchmark cases 5–7.
No foamForNuclear dependency — imported by case.py.
Case 8 (temperature cycling) lives in ../pyc_cycling/.
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

# Polynomial swelling coefficients for PyCCorrelation (IAEA TECDOC-1674 Table 9.x).
PYC_SWELLING_COEFFICIENTS = {
    "A": {   # cases 5 and 6
        "radial":     [-2.22642e-2,  2.00861e-2, -7.77024e-3,  1.36334e-3, 0.0, 0.0],
        "tangential": [-1.91253e-2,  2.63307e-3,  1.69251e-3, -3.53804e-4, 0.0, 0.0],
    },
    "B": {   # case 7
        "radial":     [-2.12522e-2,  1.83715e-2, -5.05553e-3,  7.27026e-4, 0.0, 0.0],
        "tangential": [-1.79113e-2, -3.42182e-3,  5.03465e-3, -8.88086e-4, 0.0, 0.0],
    },
}

# Fast flux: constant 3.4722e13 n/cm²/s over 1000 days
PYC_FAST_FLUX = dict(
    timePoints=[0.0, 8.64e7],
    fastFlux=[3.4722e13, 3.4722e13],
    materials=["IPyC", "SiC", "OPyC"],
)

# Layer thicknesses in μm (same for all cases): IPyC=40, SiC=35, OPyC=40
_LAYERS = [("IPyC", 40.0, 8), ("SiC", 35.0, 8), ("OPyC", 40.0, 8)]
_MAT_NAMES = ["IPyC", "SiC", "OPyC"]
_T_UNIFORM = dict(T_int=1273.0, T_inner=1273.0, T_outer=1273.0)
_TIME_PYC  = dict(endTime=8.64e7, deltaT=8.64e3, minDeltaT=1.0, maxDeltaT=4.32e5)
_CREEP_A   = dict(pyc_Tref=298.15, pyc_swelling_set="A", pyc_creep_coefficient=2.71e-4)
_CREEP_B   = dict(pyc_Tref=298.15, pyc_swelling_set="B", pyc_creep_coefficient=2.71e-4)

CASES = {
    # case_5: pressure 0 → 15.54 MPa, swelling set A
    "case_5": dict(
        inner_radius=275.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        **_T_UNIFORM,
        p_inner_list=[[0.0, 0.0], [8.64e7, 15.54e6]], p_outer=0.1e6,
        **_TIME_PYC, **_CREEP_A, restart_from=None,
    ),
    # case_6: pressure 0 → 26.2 MPa, swelling set A
    "case_6": dict(
        inner_radius=350.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        **_T_UNIFORM,
        p_inner_list=[[0.0, 0.0], [8.64e7, 26.2e6]], p_outer=0.1e6,
        **_TIME_PYC, **_CREEP_A, restart_from=None,
    ),
    # case_7: same as case_6 but swelling set B
    "case_7": dict(
        inner_radius=350.0, layers=_LAYERS, mat_names=_MAT_NAMES,
        **_T_UNIFORM,
        p_inner_list=[[0.0, 0.0], [8.64e7, 26.2e6]], p_outer=0.1e6,
        **_TIME_PYC, **_CREEP_B, restart_from=None,
    ),
}
