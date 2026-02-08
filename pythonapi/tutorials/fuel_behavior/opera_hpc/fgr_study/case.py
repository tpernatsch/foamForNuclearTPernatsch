# case.py
from foamForNuclear import (
    mesh, fields, boundaryConditions as bc, offbeat_lib as offbeat, functions)
from foamForNuclear.offbeat_lib import fast_flux, materials, behaviour, laws
from foamForNuclear.case import OffbeatCase
from foamForNuclear.offbeat_lib.materials.models import ActinideDict
from foamForNuclear.common import SciantixDict

# Geometry and mesh input
mm = 1e-3

fuel_length = 3710 * mm
plenum_length = 323.5 * mm

fuel_ri     = 0.0 * mm
fuel_ro     = 4.35 * mm

clad_ri     = 4.436 * mm
clad_ro     = 5.036 * mm

fuel_nr = 30
clad_nr = 10
slices  = 1

# default values (just for single-case runs)
h_gap = 2500    # W/m2K
T0 = 293        # K
T_cold = 573        # K
time_pts_hrs = [0, 40e3]
fast_flux = 1e13

coolant_pressure_list = [
    (0, 7e6),
    (365, 7e6)]

def build_case(
    case_name: str,
    lhgr: float,
    r_grain:float,
    ):
    # %% 1) Mesh
    rod_mesh = mesh.rod_1d(
        fuel_length=fuel_length,
        plenum_length=plenum_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        fuel_nr=fuel_nr,
        clad_nr=clad_nr,
        fuel_grading_nr=1,
        slices=slices,
    )

    # %% 2) Fields
    T = fields.Temperature(
        internalField=T0,
        boundaryField=dict(
            fuelOuter = bc.FuelRodGap(value=T0, roughness=2e-6),
            cladInner = bc.FuelRodGap(value=T0, roughness=0.5e-6),
            cladOuter = bc.FixedValue(value=T_cold),
            )
        )

    D = fields.Displacement(
        internalField=[0, 0, 0],
        boundaryField=dict(
            fuelOuter = bc.GapContact(value=[0, 0, 0], penaltyFactor=0.1, relaxInterfacePressure=0.1),
            cladInner = bc.GapContact(value=[0, 0, 0], penaltyFactor=0.1, relaxInterfacePressure=0.1),
            cladOuter = bc.CoolantPressure(value=[0, 0, 0], coolantPressureList=coolant_pressure_list),
            )
        )
    
    gap_gas = fields.GapGas(gasPressure=2.0e6)

    # %% 3) Materials       
    uo2 = materials.UO2(
        name="fuel", 
        densityFraction=0.95, 
        theoreticalDensity=10960, 
        rGrain=r_grain,
        dishFraction=0.02
        )
    uo2.isotopes = ({
        "U": ActinideDict(
            ratioOverMetal = 1.0,
            massNumbers = [235, 238],
            weightFractions = [0.05, 0.95]
    )})
    uo2.damage = materials.damage.IsotropicCracking()
    
    uo2.densification = behaviour.densification.Empirical(
        densificationDensityChange=1,
        densificationTimeConstant=1000
    )
    uo2.relocation = behaviour.relocation.UO2Frapcon(
        GapCold=(clad_ri-fuel_ro)*2,
        DiamCold=fuel_ro*2,
        outerPatch="fuelOuter"
    )
    uo2.swelling = behaviour.swelling.UO2Frapcon()

    zircaloy = materials.Zircaloy(name="cladding")
    zircaloy.constitutiveLaw = laws.MisesPlasticCreep(
        yieldStress=laws.yield_stress.Constant(sigmaY=1e9),
        creep=laws.creep.Limback(relax=0.9)
    )

    # %% 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=rod_mesh, caseFolder=case_name)
    
    case.thermalSolver    = offbeat.thermal_solver.SolidConduction()
    case.mechanicsSolver    = offbeat.mechanics_solver.SmallStrain()
    case.rheology = offbeat.rheology.ModifiedPlaneStrain(
        springModulus=3500,
        coolantPressureList=coolant_pressure_list
    )
    case.gapGasModel = offbeat.gap_gas.Frapcon(
        gapPatches=["fuelOuter", "cladInner"],
        topFuelPatches=["fuelTop"],
        bottomFuelPatches=["fuelBottom"],
    )
    
    case.burnup = offbeat.burnup.FromPower()

    case.heatSource = offbeat.heat_source.TimeDependentLhgr(
        timePoints=time_pts_hrs,
        lhgr=[lhgr, lhgr],
        materials=("fuel",),
    ) 
    case.fastFlux = offbeat.fast_flux.TimeDependentAxialProfile(
        timePoints=time_pts_hrs,
        fastFlux=[fast_flux, fast_flux],
        materials=("fuel","cladding")
    ) 

    case.fissionGasRelease = offbeat.fgr.SCIANTIX(
        SCIANTIX=SciantixDict(
            iverification = 0,
            igrain_growth = 1,
            iinert_gas_behavior = 1,
            igas_diffusion_coefficient = 2,
            iintra_bubble_evolution = 1,
            ibubble_radius = 1,
            iresolution_rate = 1,
            itrapping_rate = 1,
            inucleation_rate = 1,
            isolver = 1,
            iformat_output = 1,
            igrain_boundary_vacancy_diffusion_coefficient = 2,
            igrain_boundary_behaviour = 1,
            igrain_boundary_micro_cracking = 1,
            igrain_recrystallization = 0,
            ifuel_reactor_type = 0,
            igas_effective_coefficient = 0,
            igas_sweeping = 0,
        )
    )
    
    case.materials = [uo2, zircaloy]
    
    case.fields = [T, gap_gas, D]    

    # %% 5) Settings + functions
    case.stressAnalysis.maxOuterIter = 100

    case.settings.endTime       = 40000
    case.settings.deltaT        = 0.1
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 100
    case.settings.userTime      = "hours"

    case.settings.adjustableTimeStep         = True
    case.settings.maxDeltaT                  = 100
    case.settings.minDeltaT                  = 0.01
    case.settings.maxRelativeDeltaTIncrease  = 1e9
    case.settings.minRelativeDeltaTDecrease  = 1e9
    case.settings.maxRelativePowerIncrease   = 1e9
    case.settings.maxRelativePowerDecrease   = 1e9

    # %% 6) FunctionObjects
    fgr = functions.FGR(
        name="fgr",
    )

    case.functions = [fgr,]
    
    return case, rod_mesh


if __name__ == "__main__":
    case, rod_mesh = build_case()