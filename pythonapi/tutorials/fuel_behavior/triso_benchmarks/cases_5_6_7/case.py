from foamForNuclear import mesh, fields, functions, boundaryConditions as bc, offbeat_lib as offbeat
from foamForNuclear.case import OffbeatCase
from foamForNuclear.offbeat_lib import materials, behaviour, laws
from foamForNuclear.offbeat_lib.materials import properties as mat_props
from foamForNuclear.timeProfile import OffbeatTimeProfile
from foamForNuclear.common import Vector

from _config import MAT as _MAT, CASES as _CASES, PYC_SWELLING_COEFFICIENTS, PYC_FAST_FLUX

um = 1e-6


def _thicknesses_to_outer_radii(inner_radius, layers):
    r = inner_radius
    result = []
    for name, thickness, n in layers:
        r += thickness
        result.append((name, r, n))
    return result


def _build_pyc_material(name, Tref, swelling_set, creep_coefficient):
    m  = _MAT[name]
    sw = PYC_SWELLING_COEFFICIENTS[swelling_set]
    return materials.PyC(
        name=name, Tref=Tref,
        density=mat_props.density.Constant(value=m["density"]),
        conductivity=mat_props.conductivity.Constant(value=m["conductivity"]),
        heatCapacity=mat_props.heat_capacity.Constant(value=m["heatCapacity"]),
        emissivity=mat_props.emissivity.Constant(value=m["emissivity"]),
        YoungModulus=mat_props.young_modulus.Constant(value=m["YoungModulus"]),
        PoissonRatio=mat_props.poisson_ratio.Constant(value=m["PoissonRatio"]),
        thermalExpansion=mat_props.thermal_expansion.Constant(
            value=m["thermalExpansion"], Tref=Tref,
        ),
        swelling=behaviour.swelling.PyCCorrelation(
            radialCoefficients=sw["radial"],
            tangentialCoefficients=sw["tangential"],
        ),
        constitutiveLaw=laws.MisesPlasticCreep(
            yieldStress=laws.yield_stress.Constant(sigmaY=1e60),
            creep=laws.creep.ConstantPrincipalStress(creepCoefficient=creep_coefficient),
        ),
    )


def build_case(case_id: str = "case_5"):
    if case_id not in _CASES:
        raise ValueError(f"Unknown case '{case_id}'. Available: {list(_CASES.keys())}")

    cfg = _CASES[case_id]

    sphere_mesh = mesh.sphere_1d(
        wedge_angle=0.25, scale=um,
        inner_radius=cfg["inner_radius"],
        layers=_thicknesses_to_outer_radii(cfg["inner_radius"], cfg["layers"]),
    )

    case = OffbeatCase(case_id, mesh=sphere_mesh)

    T = fields.Temperature(
        internalField=cfg["T_int"],
        boundaryField=dict(
            inner=bc.FixedValue(value=cfg["T_inner"]),
            outer=bc.FixedValue(value=cfg["T_outer"]),
        )
    )

    D = fields.Displacement(
        internalField=[0, 0, 0],
        boundaryField=dict(
            inner=bc.TractionDisplacement(
                value=[0, 0, 0], traction=Vector(0, 0, 0),
                pressureList=OffbeatTimeProfile(type="table", values=cfg["p_inner_list"]),
            ),
            outer=bc.TractionDisplacement(
                value=[0, 0, 0], traction=Vector(0, 0, 0), pressure=cfg["p_outer"],
            ),
        )
    )

    case.fields += [T, D]

    Tref   = cfg["pyc_Tref"]
    sw_set = cfg["pyc_swelling_set"]
    creep  = cfg["pyc_creep_coefficient"]

    case.materials = [
        _build_pyc_material(name, Tref, sw_set, creep) if name in ("IPyC", "OPyC")
        else materials.Constant(name=name, **_MAT[name])
        for name in cfg["mat_names"]
    ]

    case.thermalSolver   = offbeat.thermal_solver.SolidConduction()
    case.mechanicsSolver = offbeat.mechanics_solver.SmallStrain(cylindricalStress=True)
    case.rheology        = offbeat.rheology.Standard(thermalExpansion=False)

    case.fastFlux = offbeat.fast_flux.TimeDependentAxialProfile(
        timePoints=PYC_FAST_FLUX["timePoints"],
        fastFlux=PYC_FAST_FLUX["fastFlux"],
        materials=PYC_FAST_FLUX["materials"],
    )

    r_in  = cfg["inner_radius"] * um
    r_out = r_in + sum(t for _, t, _ in cfg["layers"]) * um

    case.functions = [functions.Graph(
        name="radialProfile", start=[r_in - 1e-6, 0.0, 0.0], end=[r_out + 1e-6, 0.0, 0.0],
        graph_type="midPoint",
        fields=["T", "sigma"], nPoints=50, axis="x",
        writeControl="timeStep",
        interpolationScheme="cellPointFace",
    )]

    case.settings.endTime       = cfg["endTime"]
    case.settings.deltaT        = cfg["deltaT"]
    case.settings.writeControl  = "runTime"
    case.settings.writeInterval = 4.32e6

    case.settings.adjustableTimeStep        = True
    case.settings.minDeltaT                 = cfg["minDeltaT"]
    case.settings.maxDeltaT                 = cfg["maxDeltaT"]
    case.settings.maxRelativeDeltaTIncrease = 1e9
    case.settings.minRelativeDeltaTDecrease = 1e9
    case.settings.maxBurnupIncrease         = 0.1
    case.settings.maxAverageCreep           = 1e-4
    case.settings.maxMaximumCreep           = 1e-4

    return case, sphere_mesh


if __name__ == "__main__":
    import sys
    build_case(sys.argv[1] if len(sys.argv) > 1 else "case_5")
