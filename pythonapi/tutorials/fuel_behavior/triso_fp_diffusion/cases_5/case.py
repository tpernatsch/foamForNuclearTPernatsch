from foamForNuclear import mesh, fields, boundaryConditions as bc, offbeat_lib as offbeat
from foamForNuclear.case import OffbeatCase
from foamForNuclear.common import Table
from foamForNuclear.offbeat_lib import materials
from foamForNuclear.offbeat_lib.materials.models import ActinideDict
from foamForNuclear.offbeat_lib import element_transport as et
from foamForNuclear.offbeat_lib.element_transport import diffCoef_models as dm
from foamForNuclear.offbeat_lib.element_transport import yield_models as ym
from foamForNuclear.functions import FpRelease

from _config import (MAT, CASES,
                     N_CYCLES, CYCLE_DURATION,
                     T_IRRAD_LOW, T_IRRAD_HIGH, T_HEATING,
                     VHGR, FAST_FLUX_RATE)

um = 1e-6


def _cyclic_T_table(cfg: dict) -> list:
    """Cyclic irradiation: linear ramp T_LOW→T_HIGH each cycle, instant drop between cycles."""
    rows = [[0, T_IRRAD_LOW]]
    for k in range(N_CYCLES):
        t_end = (k + 1) * CYCLE_DURATION
        rows.append([t_end, T_IRRAD_HIGH])
        if k < N_CYCLES - 1:
            rows.append([t_end + 3600, T_IRRAD_LOW])
    if cfg.get("t_heating") is not None:
        t_h = cfg["t_heating"]
        rows.append([t_h + 3600, T_HEATING])
        rows.append([cfg["endTime"], T_HEATING])
    return rows


def build_case(case_id: str = "case_5b"):
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
    outer_bc = bc.UniformFixedValue(uniformValue=Table(_cyclic_T_table(cfg)))

    T = fields.Temperature(
        internalField=T_IRRAD_LOW,
        boundaryField=dict(outer=outer_bc),
    )

    # Start from zero: Cs and Ag are produced during irradiation via yield model
    N_Cs = fields.Field(
        name="N_Cs",
        dimensions=fields.Dimension(moles=1, length=-3),
        internalField=0.0,
        boundaryField=dict(outer=bc.FixedValue(value=0.0)),
    )

    N_Ag = fields.Field(
        name="N_Ag",
        dimensions=fields.Dimension(moles=1, length=-3),
        internalField=0.0,
        boundaryField=dict(outer=bc.FixedValue(value=0.0)),
    )

    case.fields += [T, N_Cs, N_Ag]

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
    case.thermalSolver = offbeat.thermal_solver.SolidConduction()
    case.rheology      = offbeat.rheology.Standard(thermalExpansion=False)

    t_end = cfg["endTime"]
    t_heat_start = cfg.get("t_heating")
    zone_names = [layer[0] for layer in cfg["layers"]]

    if t_heat_start is not None:
        # case_5b: fission active during irradiation, off during heating
        time_pts  = [0, t_heat_start, t_heat_start + 3600, t_end]
        vhgr_vals = [VHGR, VHGR, 0.0, 0.0]
        flux_vals = [FAST_FLUX_RATE, FAST_FLUX_RATE, 0.0, 0.0]
    else:
        # case_5a: fission active throughout
        time_pts  = [0, t_end]
        vhgr_vals = [VHGR, VHGR]
        flux_vals = [FAST_FLUX_RATE, FAST_FLUX_RATE]

    case.heatSource = offbeat.heat_source.TimeDependentVhgr(
        timePoints=time_pts,
        vhgr=vhgr_vals,
        materials=["Kernel"],
    )
    case.burnup = offbeat.burnup.FromPower()

    case.fastFlux = offbeat.fast_flux.TimeDependentAxialProfile(
        timePoints=time_pts,
        fastFlux=flux_vals,
        materials=zone_names,
    )

    case.elementTransportSolver = et.FpDiffusion(
        fissionProducts=["Cs", "Ag"],
        diffusion={
            "Kernel":    dm.UO2Arrhenius(),
            "Buffer":    dm.SingleArrhenius(d1=1e-8, q1=0.0),
            "IPyC|OPyC": dm.PyCArrhenius(),
            "SiC":       dm.SiCArrhenius(),
        },
        yield_={
            "Kernel": {
                "Cs": ym.EnrichmentDependent(yieldLEU=0.14, yieldHEU=0.16),
                "Ag": ym.BurnupDependent(yieldCoeff=1.31625e-3, burnupExponent=0.55734),
            },
        },
    )

    # ── Function objects ──────────────────────────────────────────────────────
    case.add_function_object(FpRelease(name="writeFpRelease"))

    # ── Settings ──────────────────────────────────────────────────────────────
    case.settings.endTime       = t_end
    case.settings.deltaT        = cfg["deltaT"]
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 480
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
    build_case(sys.argv[1] if len(sys.argv) > 1 else "case_5b")
