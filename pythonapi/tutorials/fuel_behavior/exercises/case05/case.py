# %% Imports

from foamForNuclear import (
    mesh, functions, fields, boundaryConditions as bc, offbeat_lib as offbeat)
from foamForNuclear.case import OffbeatCase
from foamForNuclear.offbeat_lib import materials, behaviour, laws
from foamForNuclear.offbeat_lib.materials.models import ActinideDict


# %% Parameters (edit here)
mm = 1e-3

fuel_length      = 110.0 * mm
plenum_length    = 10.0  * mm

fuel_ri          = 0.0   * mm
fuel_ro          = 4.5   * mm

gap              = 0.06  * mm
clad_thickness   = 0.7   * mm

clad_ri          = fuel_ro + gap
clad_ro          = clad_ri + clad_thickness

fuel_nr = 30
clad_nr = 10
slices  = 1

T0    = 293.0     # K
Tcold = 573.0     # K

# Operating conditions
time_pts_days = [0, 0.04, 364, 364.04, 365]
lhgr_W_per_m  = [0, 40000, 40000, 0, 0]

fast_flux     = [0, 1e13, 1e13, 0, 0]
    
coolant_pressure_list = [
    (0, 1e7),
    (365, 1e7)]

def build_case():
    # %% 1) Geometry and mesh
    rod_mesh = mesh.rod_1d(
        fuel_length=fuel_length,
        plenum_length=plenum_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        fuel_nr=fuel_nr,
        clad_nr=clad_nr,
        slices=slices,
    )

    # %% 2) Field & BCs
    T = fields.Temperature(
        internalField=T0,
        boundaryField=dict(
            fuelOuter = bc.FuelRodGap(value=T0, roughness=2e-6),
            cladInner = bc.FuelRodGap(value=T0, roughness=0.5e-6),
            cladOuter = bc.FixedValue(value=Tcold),
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

    # %% 3) Materials (constant k)
    uo2 = materials.UO2(
        name="fuel",
        rGrain=2.8e-5)
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
        GapCold=0.12e-3,
        DiamCold=9e-3,
        outerPatch="fuelOuter"
    )

    zircaloy = materials.Zircaloy(name="cladding")
    zircaloy.constitutiveLaw = laws.MisesPlasticCreep(
        yieldStress=laws.yield_stress.Constant(sigmaY=400e6),
        creep=laws.creep.Limback(relax=0.9)
    )

    # %% 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=rod_mesh, caseFolder="case")

    case.burnup = offbeat.burnup.FromPower()

    case.thermalSolver = offbeat.thermal_solver.SolidConduction()
    case.mechanicsSolver = offbeat.mechanics_solver.SmallStrain()
    case.rheology = offbeat.rheology.ModifiedPlaneStrain(
        springModulus=3500,
        coolantPressureList=coolant_pressure_list
    )
    case.gapGasModel = offbeat.gap_gas.Frapcon(
        gapPatches=["fuelOuter", "cladInner"],
        topFuelPatches=["fuelTop"],
        bottomFuelPatches=["fuelBottom"],
    )

    case.heatSource = offbeat.heat_source.TimeDependentLhgr(
        timePoints=time_pts_days,
        lhgr=lhgr_W_per_m,
        materials=("fuel",)
    )
    case.fastFlux = offbeat.fast_flux.TimeDependentAxialProfile(
        timePoints=time_pts_days,
        fastFlux=fast_flux,
        materials=("fuel","cladding")
    )

    case.materials = [uo2, zircaloy]

    case.fields = [T, gap_gas, D]    

    # %%  5) Settings + functions
    case.stressAnalysis.maxOuterIter = 100

    case.settings.endTime       = 365
    case.settings.deltaT        = 0.01
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 1
    case.settings.userTime      = "days"

    case.settings.adjustableTimeStep         = True
    case.settings.maxDeltaT                  = 7
    case.settings.minDeltaT                  = 0.01
    case.settings.maxRelativeDeltaTIncrease  = 1e9
    case.settings.minRelativeDeltaTDecrease  = 1e9
    case.settings.maxRelativePowerIncrease   = 1e9
    case.settings.maxRelativePowerDecrease   = 1e9

    # %% 6) FunctionObjects
    fuelCenterline = functions.Probes(
        name="fuelCenterline",
        fields=['T', 'Bu'],
        probeLocations=[[fuel_ri, 0.0, fuel_length/2]],
        enabled=True
    )

    fuelOuter = functions.PatchProbes(
        name="fuelOuter",
        fields=['T', 'Bu', 'gapWidth'],
        probeLocations=[[fuel_ri, 0.0, fuel_length/2]],
        patches=["fuelOuter",],
        enabled=True
    )

    radialProfile = functions.Graph(
        name="radialProfile",
        start=[fuel_ri, 0.0, fuel_length/2],
        end=[fuel_ro, 0.0,  fuel_length/2],
        nPoints=30,
        fields=['T', 'sigma'],
        interpolationScheme="cell"
    )

    case.functions = [fuelCenterline, radialProfile, fuelOuter]

    return case, rod_mesh

if __name__ == "__main__":
    case, rod_mesh = build_case()