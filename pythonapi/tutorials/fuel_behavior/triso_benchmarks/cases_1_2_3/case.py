from foamForNuclear import mesh, fields, functions, boundaryConditions as bc, offbeat_lib as offbeat
from foamForNuclear.case import OffbeatCase
from foamForNuclear.offbeat_lib import materials
from foamForNuclear.common import Vector

from _config import MAT as _MAT, CASES as _CASES

um = 1e-6


def build_case(case_id: str = "case_1"):
    if case_id not in _CASES:
        raise ValueError(f"Unknown case '{case_id}'. Available: {list(_CASES.keys())}")

    cfg = _CASES[case_id]

    sphere_mesh = mesh.sphere_1d(
        wedge_angle=0.25, scale=um,
        inner_radius=cfg["inner_radius"], layers=cfg["layers"],
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
                value=[0, 0, 0], traction=Vector(0, 0, 0), pressure=cfg["p_inner"],
            ),
            outer=bc.TractionDisplacement(
                value=[0, 0, 0], traction=Vector(0, 0, 0), pressure=cfg["p_outer"],
            ),
        )
    )

    case.fields += [T, D]
    case.materials = [materials.Constant(name=name, **_MAT[name]) for name in cfg["mat_names"]]

    case.thermalSolver   = offbeat.thermal_solver.SolidConduction()
    case.mechanicsSolver = offbeat.mechanics_solver.SmallStrain(cylindricalStress=True)
    case.rheology        = offbeat.rheology.Standard(thermalExpansion=False)

    r_in  = cfg["inner_radius"] * um
    r_out = cfg["layers"][-1][1] * um

    case.functions = [functions.Graph(
        name="radialProfile", start=[r_in - 1e-6, 0.0, 0.0], end=[r_out + 1e-6, 0.0, 0.0],
        graph_type="midPoint",
        fields=["T", "sigma"], nPoints=50, axis="x",
        writeControl="timeStep", interpolationScheme="cellPointFace",
    )]

    case.settings.endTime       = cfg["endTime"]
    case.settings.deltaT        = cfg["deltaT"]
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 2

    case.settings.adjustableTimeStep        = True
    case.settings.minDeltaT                 = cfg["minDeltaT"]
    case.settings.maxRelativeDeltaTIncrease = 1e9
    case.settings.minRelativeDeltaTDecrease = 1e9
    case.settings.maxBurnupIncrease         = 0.1
    case.settings.maxAverageCreep           = 1e-4
    case.settings.maxMaximumCreep           = 1e-4

    return case, sphere_mesh


if __name__ == "__main__":
    import sys
    build_case(sys.argv[1] if len(sys.argv) > 1 else "case_1")
