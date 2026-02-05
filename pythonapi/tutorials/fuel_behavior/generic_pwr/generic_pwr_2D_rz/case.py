# %% Imports

from foamForNuclear import (
    mesh, functions, fields, boundaryConditions as bc, offbeat_lib as offbeat)
from foamForNuclear.case import OffbeatCase
from foamForNuclear.offbeat_lib import materials, behaviour, laws
from foamForNuclear.offbeat_lib.materials.models import ActinideDict
from foamForNuclear.offbeat_lib import profiles
from foamForNuclear.common import SciantixDict


# %% Parameters (edit here)
mm = 1e-3

# We only represent a portion of the full rod geometry in this example
fuel_length      = 300.0 * mm
plenum_length    = 20.0  * mm

fuel_ri          = 0.0   * mm
fuel_ro          = 4.5   * mm

clad_ri          = 4.565 * mm
clad_ro          = 5.315 * mm

fuel_nr = 30
fuel_nz = 30

clad_nr = 10
clad_nz = 32

# Operating conditions
T0    = 300.0     # K
time_pts_hrs  = [0, 1, 35e3]
lhgr_W_per_m  = [0, 20e3, 10e3]
fast_flux     = [0, 2e13, 1e13]
    
coolant_pressure_list = [
    (0, 1e5),
    (1, 1e7),
    (35e3, 1e7)]

cladding_temperature_list = [
    (0, 300),
    (1, 600),
    (35e3, 600)]

# Lhgr axial profile
axial_profile_time_points = [
    0.0000E+00, 4.3889E+03, 8.7500E+03, 8.7778E+03, 1.3139E+04, 1.7528E+04, 1.7556E+04, 2.1889E+04, 2.6278E+04, 2.6306E+04, 3.0556E+04, 3.5000E+04
    ]

axial_profile_locations = [
    0, 0.090909090909091, 0.181818181818182, 0.272727272727273,
    0.363636363636364, 0.454545454545455, 0.545454545454546,
    0.636363636363636, 0.727272727272727, 0.818181818181818,
    0.909090909090909, 1]

axial_profile_data = [
    [0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ],
    [ 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ],
    [ 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ],
    
    [0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ],
    [ 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ],
    [ 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ],    
    
    [0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ],
    [ 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ],
    [ 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ],    
    
    [0.469787352107652, 0.913850849204304, 1.19921124080937, 1.32281654412493, 1.3457064151093, 1.29840068174161, 1.21141917200103, 1.08628787728651, 0.942844685784499, 0.767355674904378, 0.541508947858657, 0.271408470243168 ],
    [ 0.625964954658454, 1.03050615457391, 1.13900379439798, 1.13280392926518, 1.10335456988436, 1.08630494076915, 1.07235524422034, 1.06305544652113, 1.04135591855632, 0.977807300945077, 0.811960908642572, 0.457018629789547 ],
    [ 0.807036867064403, 1.11329355973688, 1.10701137116924, 1.05518331548621, 1.01749018408036, 1.00178471266126, 0.997073071235531, 1.00335525980317, 1.02063127836418, 1.023772372648, 0.942103921268676, 0.629565040028555 ],]

def build_case():
    # %% 1) Geometry and mesh
    rod_mesh = mesh.rod_2d_rz(
        fuel_length=fuel_length,
        plenum_length=plenum_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        fuel_nr=fuel_nr,
        fuel_nz=fuel_nz,
        clad_nr=clad_nr,
        clad_nz=clad_nz,
        updateAMI=False
    )
    
    # This test case represents only a portion of the full rod (a rodlet), 
    # therefore we tranlsate it to mid rod position
    rod_mesh.transforms.append(("translate", (0.0, 0.0, 1.35)))

    # %% 2) Field & BCs
    T = fields.Temperature(
        internalField=T0,
        boundaryField=dict(
            fuelOuter = bc.FuelRodGap(value=T0, roughness=2e-6, relax=0.1),
            cladInner = bc.FuelRodGap(value=T0, roughness=0.5e-6, relax=0.1),
            cladOuter = bc.UniformFixedValue(uniformValue=cladding_temperature_list)
            )
        )

    D = fields.Displacement(
        internalField=[0, 0, 0],
        boundaryField=dict(
            fuelOuter = bc.GapContact(value=[0, 0, 0], penaltyFactor=0.1, relaxInterfacePressure=0.1),
            cladInner = bc.GapContact(value=[0, 0, 0], penaltyFactor=0.1, relaxInterfacePressure=0.1),
            cladOuter = bc.CoolantPressure(value=[0, 0, 0], coolantPressureList=coolant_pressure_list),
            fuelBottom = bc.FixedDisplacementZeroShear(value=[0, 0, 0]),
            cladBottom = bc.FixedDisplacementZeroShear(value=[0, 0, 0]),
            fuelTop    = bc.PlenumSpringPressure(
                value=[0, 0, 0],
                fuelTopPatches=["fuelTop"],
                topCapInnerPatches=["cladTop"],
                springModulus=3.5e3,
                planeStrain=True,
            ),
            cladTop    = bc.TopCladRingPressure(
                value=[0, 0, 0],
                topCapOuterRadius=5.325e-3,
                topCapInnerRadius=4.575e-3,
                coolantPressureList=coolant_pressure_list,
                planeStrain=True,
            )
        ))

    gap_gas = fields.GapGas(gasPressure=2.25e6)
    
    neutron_flux = fields.NeutronFlux0(
        internalField=0.0,
        boundaryField=dict(
            fuelOuter = bc.FixedValue(value=1),
        ),
    )

    # %% 3) Materials (constant k)
    uo2 = materials.UO2(name="fuel", rGrain=5e-6)
    uo2.isotopes = ({
        "U": ActinideDict(
            ratioOverMetal = 1.0,
            massNumbers = [235, 238],
            weightFractions = [0.045, 1 - 0.045]
    )})
    uo2.damage = materials.damage.IsotropicCracking()
    
    uo2.densification = behaviour.densification.UO2Frapcon(
        resinteringDensityChange=1
    )
    uo2.relocation = behaviour.relocation.UO2Frapcon(
        GapCold=2*(clad_ri-fuel_ro),
        DiamCold=2*fuel_ro,
        outerPatch="fuelOuter",
        recoveryFraction=0.5
    )
    uo2.swelling = behaviour.swelling.UO2Frapcon()

    zircaloy = materials.Zircaloy(name="cladding")
    zircaloy.constitutiveLaw = laws.MisesPlasticCreep(
        yieldStress=laws.yield_stress.Constant(sigmaY=1e9),
        creep=laws.creep.Limback(relax=0.9)
    )

    # %% 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=rod_mesh, caseFolder="case")

    case.neutronicsSolver = offbeat.neutronics_solver.Diffusion()
    case.burnup = offbeat.burnup.Lassmann(convergencePrecision=1e-6)

    case.thermalSolver = offbeat.thermal_solver.SolidConduction()
    case.mechanicsSolver = offbeat.mechanics_solver.SmallStrain()
    
    case.gapGasModel = offbeat.gap_gas.Frapcon(
        gapPatches=["fuelOuter", "cladInner"],
        topFuelPatches=["fuelTop"],
        bottomFuelPatches=["fuelBottom"],
    )

    # NOTE: one can still use the full rod axial profile by using the 
    # useExternalReferenceDimensions keyword
    case.heatSource = offbeat.heat_source.TimeDependentLhgr(
        timePoints=time_pts_hrs,
        lhgr=lhgr_W_per_m,
        materials=("fuel",),
        axialProfile=profiles.axial_profile.TimeDependentTabulated(
            timePoints=axial_profile_time_points,
            axialLocations=axial_profile_locations,
            data=axial_profile_data
        ),
        radialProfile=profiles.radial_profile.FromBurnup(),
        useExternalReferenceDimensions=True,
        zMin=0.0,
        zMax=3.0,
        volumeFraction=0.1
    )

    case.fastFlux = offbeat.fast_flux.TimeDependentAxialProfile(
        timePoints=time_pts_hrs,
        fastFlux=fast_flux,
        materials=("fuel","cladding")
    )

    case.fissionGasRelease = offbeat.fgr.SCIANTIX(
        SCIANTIX=SciantixDict(
            iverification = 0,
            igrain_growth = 1,
            iinert_gas_behavior = 1,
            igas_diffusion_coefficient = 1,
            iintra_bubble_evolution = 1,
            ibubble_radius = 1,
            iresolution_rate = 1,
            itrapping_rate = 1,
            inucleation_rate = 1,
            isolver = 1,
            iformat_output = 1,
            igrain_boundary_vacancy_diffusion_coefficient = 1,
            igrain_boundary_behaviour = 1,
            igrain_boundary_micro_cracking = 1,
            igrain_recrystallization = 0,
            ifuel_reactor_type = 0,
            igas_effective_coefficient = 1,
            igas_sweeping = 1,
        )
    )

    case.materials = [uo2, zircaloy]

    case.fields = [T, gap_gas, D, neutron_flux]    

    # %%  5) Settings + functions
    case.solver.fvSolution.solvers["D"].solver = "PBiCGStab"
    case.stressAnalysis.nCorrectors = dict(
        default = 1,
        D = 3
    )
    case.stressAnalysis.maxOuterIter = 100

    case.settings.endTime       = 35e3
    case.settings.deltaT        = 0.001
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
    case.settings.maxBurnupIncrease          = 0.1
    case.settings.maxAverageCreep            = 1e-4
    case.settings.maxMaximumCreep            = 1e-4

    # %% 6) FunctionObjects
    offset=1.35
    fuelCenterline = functions.Probes(
        name="fuelCenterline",
        fields=['T', 'Bu'],
        probeLocations=[[fuel_ri, 0.0, fuel_length/2+offset]],
        enabled=True
    )

    fuelOuter = functions.PatchProbes(
        name="fuelOuter",
        fields=['T', 'Bu', 'gapWidth', 'interfaceP', 'hGap'],
        probeLocations=[[fuel_ri, 0.0, fuel_length/2+offset]],
        patches=["fuelOuter",],
        enabled=True
    )

    dz = fuel_length/fuel_nz
    slice_id = 4
    radialProfile = functions.Graph(
        name="radialProfile",
        start=[fuel_ri, 0.0, (slice_id + 0.5)*dz+offset ],
        end=[fuel_ro, 0.0,  (slice_id + 0.5)*dz+offset],
        fields=['T', 'sigma'],
        nPoints=100,
        writeControl="timeStep",
        interpolationScheme="cellPatchConstrained"
    )

    fgr = functions.FGR(
        name="fgr",
    )

    case.functions = [fuelCenterline, radialProfile, fuelOuter, fgr]

    return case, rod_mesh

if __name__ == "__main__":
    case, rod_mesh = build_case()