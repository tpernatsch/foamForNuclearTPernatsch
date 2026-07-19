from foamForNuclear import mesh, fields, boundaryConditions as bc, offbeat_lib as offbeat
from foamForNuclear.case import OffbeatCase
from foamForNuclear.common import Table
from foamForNuclear.offbeat_lib import materials
from foamForNuclear.offbeat_lib.materials.models import ActinideDict
from foamForNuclear.offbeat_lib import element_transport as et
from foamForNuclear.offbeat_lib.element_transport import diffCoef_models as dm
from foamForNuclear.offbeat_lib.element_transport import yield_models as ym
from foamForNuclear.functions import FpRelease

from _config import MAT, CASES

um = 1e-6


def build_case(case_id: str = "case_3a"):
    if case_id not in CASES:
        raise ValueError(f"Unknown case '{case_id}'. Available: {list(CASES.keys())}")

    cfg = CASES[case_id]

    # ── Mesh ─────────────────────────────────────────────────────────────────
    sphere_mesh = mesh.sphere_1d(
        wedge_angle=0.25, scale=um,
        inner_radius=0.0, layers=cfg["layers"],
    )

    case = OffbeatCase(case_id, mesh=sphere_mesh)

    # ── Fields ───────────────────────────────────────────────────────────────
    t_switch = cfg.get("t_switch")
    if t_switch is None:
        T_init = cfg["T"]
        outer_bc = bc.FixedValue(value=cfg["T"])
    else:
        T1, T2 = cfg["T"]
        T_init = T1
        outer_bc = bc.UniformFixedValue(
            uniformValue=Table([
                [0,              T1],
                [t_switch,       T1],
                [t_switch + 1,   T2],
                [cfg["endTime"], T2],
            ])
        )

    T = fields.Temperature(
        internalField=T_init,
        boundaryField=dict(outer=outer_bc),
    )

    N_Cs = fields.Field(
        name="N_Cs",
        dimensions=fields.Dimension(moles=1, length=-3),
        internalField=0.0,
        boundaryField=dict(outer=bc.FixedValue(value=0.0)),
    )

    case.fields += [T, N_Cs]
    setfields_pairs = [(N_Cs, 1.0)]

    if cfg.get("include_kr"):
        N_Kr = fields.Field(
            name="N_Kr",
            dimensions=fields.Dimension(moles=1, length=-3),
            internalField=0.0,
            boundaryField=dict(outer=bc.FixedValue(value=0.0)),
        )
        case.fields += [N_Kr]
        setfields_pairs.append((N_Kr, 1.0))

    case.solver.setFieldsDict.add_zone_to_cell(setfields_pairs, "Kernel")

    # ── Materials ─────────────────────────────────────────────────────────────
    mat_list = []
    for name in cfg["mat_names"]:
        if name == "Kernel":
            mat_list.append(materials.UO2(
                name="Kernel",
                enrichment=0.025,
                rGrain=6e-6,
                densityFraction=0.95,
                theoreticalDensity=10960.0,
                isotopes={
                    "U": ActinideDict(
                        ratioOverMetal=1.0,
                        massNumbers=[235, 238],
                        weightFractions=[0.025, 0.975],
                    )
                },
            ))
        else:
            mat_list.append(materials.Constant(name=name, **MAT[name]))
    case.materials = mat_list

    # ── Physics ───────────────────────────────────────────────────────────────
    case.thermalSolver   = offbeat.thermal_solver.SolidConduction()
    case.rheology        = offbeat.rheology.Standard(thermalExpansion=False)

    case.heatSource = offbeat.heat_source.Constant()
    case.burnup     = offbeat.burnup.Constant()

    zone_names = [layer[0] for layer in cfg["layers"]]
    case.fastFlux = offbeat.fast_flux.TimeDependentAxialProfile(
        timePoints=[0, 1e9],
        fastFlux=[0.0, 0.0],
        materials=zone_names,
    )

    sic_model = dm.SiCArrhenius(failureTime=cfg.get("sic_failure_time"))
    pyc_model = dm.PyCArrhenius(failureTime=cfg.get("pyc_failure_time"))

    if cfg.get("include_kr"):
        case.elementTransportSolver = et.FpDiffusion(
            fissionProducts=["Cs", "Kr"],
            diffusion={
                "Kernel":    dm.UO2Arrhenius(),
                "Buffer":    dm.SingleArrhenius(d1=1e-8, q1=0.0),
                "IPyC|OPyC": pyc_model,
                "SiC":       sic_model,
            },
            yield_={
                "Kernel": {
                    "Cs": ym.EnrichmentDependent(yieldLEU=0.14, yieldHEU=0.16),
                    "Kr": ym.Constant(yield_=0.0),
                },
            },
        )
    else:
        case.elementTransportSolver = et.FpDiffusion(
            fissionProducts=["Cs"],
            diffusion={
                "Kernel":    dm.UO2Arrhenius(),
                "Buffer":    dm.SingleArrhenius(d1=1e-8, q1=0.0),
                "IPyC|OPyC": pyc_model,
                "SiC":       sic_model,
            },
            yield_={"Kernel": ym.EnrichmentDependent(yieldLEU=0.14, yieldHEU=0.16)},
        )

    # ── Function objects ──────────────────────────────────────────────────────
    case.add_function_object(FpRelease(name="writeFpRelease"))

    # ── Settings ──────────────────────────────────────────────────────────────
    case.settings.endTime       = cfg["endTime"]
    case.settings.deltaT        = cfg["deltaT"]
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 20
    case.settings.purgeWrite    = 1

    case.settings.adjustableTimeStep        = False
    case.settings.maxBurnupIncrease         = 0.1
    case.settings.maxAverageCreep           = 1e-4
    case.settings.maxMaximumCreep           = 1e-4

    # Harmonic interpolation is mandatory for multi-material diffusion.
    case.solver.fvSchemes.laplacianSchemes['default'] = 'Gauss harmonic corrected'

    return case, sphere_mesh


if __name__ == "__main__":
    import sys
    build_case(sys.argv[1] if len(sys.argv) > 1 else "case_3a")
