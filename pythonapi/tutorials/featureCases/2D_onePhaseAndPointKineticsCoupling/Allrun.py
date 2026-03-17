"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh
import foamForNuclear.porous_medium as porous


#==============================================================================*

ffn.allclean()

#==============================================================================*
# Mesh

def createMesh(region: str) -> mesh.BlockMesh:

    newMesh = mesh.BlockMesh(region=region)

    zone0 = newMesh.create_cube("zone0", 0, 0, 0, 1, 0.1, 1, nx=25, nz=25)

    inlet = mesh.Face("inlet")
    inlet.add_sub_face(zone0.bottomFace())

    outlet = mesh.Face("outlet")
    outlet.add_sub_face(zone0.topFace())

    fixedWalls = mesh.Face("fixedWalls", boundaryType="wall")
    fixedWalls.add_sub_face(zone0.leftFace())
    fixedWalls.add_sub_face(zone0.rightFace())

    frontAndBack = mesh.Face("frontAndBack", boundaryType="empty")
    frontAndBack.add_sub_face(zone0.frontFace())
    frontAndBack.add_sub_face(zone0.backFace())

    newMesh.add_boundary(inlet)
    newMesh.add_boundary(outlet)
    newMesh.add_boundary(fixedWalls)
    newMesh.add_boundary(frontAndBack)

    return(newMesh)

nMesh = createMesh(region="neutroRegion")
thMesh = createMesh(region="fluidRegion")


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.fields.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.fields.Dimension(default='neutronFlux')
defaultFlux.internalField = 1
defaultFlux.set_boundary_condition("fixedWalls", bc.FixedValue(0))
defaultFlux.set_boundary_condition("inlet", bc.FixedValue(0))
defaultFlux.set_boundary_condition("outlet", bc.FixedValue(0))
defaultFlux.set_boundary_condition("frontAndBack", bc.Empty())

T = ffn.fields.Field("T", region=thMesh.region)
T.dimensions = ffn.fields.Dimension(default='T')
T.internalField = 600
T.set_boundary_condition("fixedWalls", bc.ZeroGradient())
T.set_boundary_condition("inlet", bc.FixedValue(600))
T.set_boundary_condition("outlet", bc.ZeroGradient())
T.set_boundary_condition("frontAndBack", bc.Empty())

U = ffn.fields.Field("U", region=thMesh.region)
U.dimensions = ffn.fields.Dimension(default='U')
U.internalField = ffn.common.Vector(0, 0, 1)
U.set_boundary_condition("fixedWalls", bc.Slip())
U.set_boundary_condition("inlet", bc.FixedValue(U.internalField))
U.set_boundary_condition("outlet", bc.ZeroGradient())
U.set_boundary_condition("frontAndBack", bc.Empty())

p = ffn.fields.Field("p", region=thMesh.region)
p.dimensions = ffn.fields.Dimension(default='p')
p.internalField = 1e5
p.set_boundary_condition("fixedWalls", bc.Calculated(p.internalField))
p.set_boundary_condition("inlet", bc.Calculated(p.internalField))
p.set_boundary_condition("outlet", bc.Calculated(p.internalField))
p.set_boundary_condition("frontAndBack", bc.Empty())

p_rgh = ffn.fields.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.fields.Dimension(default='p')
p_rgh.internalField = 1e5
p_rgh.set_boundary_condition("fixedWalls", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("inlet", bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition("outlet", bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition("frontAndBack", bc.Empty())

timeFolder0.append(defaultFlux)
timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)


#==============================================================================*
# Thermal-hydraulics solver

thSolver = ffn.solvers.thermal_hydraulics.OnePhaseThermalHydraulicsSolver(
    region=thMesh.region,
    mesh=thMesh,
    removeBaffles=True,
    isSetFvSolutionToDefault=False
)

thSolver.fluid.thermophysicalProperties = ffn.thermo.SodiumConst()

thSolver.fluid.turbulenceProperties.simulationType = 'laminar'

core = porous.Structure(zones=['zone0'], volumeFraction=0.5, Dh=0.01)

core.powerModel = porous.power_models.NuclearFuelPin(
    fuelInnerRadius=0.0012,
    fuelOuterRadius=0.004715,
    cladInnerRadius=0.004865,
    cladOuterRadius=0.005365,
    fuelMeshSize=16,
    cladMeshSize=4,
    fuelRho=10480,
    fuelCp=250,
    fuelK=3,
    cladRho=7500,
    cladCp=500,
    cladK=20,
    fuelT=600,
    cladT=600,
    gapH=3000
    # gapHPowerDensity=[
    #     (   0              , 528.4782467867  ),
    #     (   795443280.75916, 3170.8694807204 )
    # ]
)

core.passiveProperties = porous.PassiveProperties(
    volumetricArea=2,
    rho=7700,
    Cp=500,
    T=600
)

thSolver.structures.append(core)

thSolver.fluid_structure.dragModels.append(
    porous.drag.ReynoldsPower(coeff=0.687, exp=-0.25, zones=['zone0'])
)
thSolver.fluid_structure.heatTransferModels.append(
    porous.heat_transfer.NusseltReynoldsPrandtlPower(
        const=4.82,
        coeff=0.0185,
        expRe=0.827,
        expPr=0.827,
        zones=['zone0']
    )
)


thSolution = ffn.numerics.fvSolution()

thSolution.add_fvSolutionSolver('p_rgh.*', ffn.numerics.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-8, relTol=0
))
thSolution.add_fvSolutionSolver('e.*', ffn.numerics.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.add_fvSolutionSolver('h.*', ffn.numerics.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.add_fvSolutionSolver('.*', ffn.numerics.fvSolutionSolver(
    solver='PBiCGStab',
    preconditioner='diagonal',
    tolerance=1e-6, relTol=0.001
))

thSolver.fvSolution = thSolution

thSolver.pimpleOptions.nOuterCorrectors = 3


thSolver.fvSchemes.divSchemes['default'] = 'none'
thSolver.fvSchemes.divSchemes['div(phi,alpha)'] = 'Gauss vanLeer'
thSolver.fvSchemes.divSchemes['div(phir,alpha)'] = 'Gauss vanLeer'
thSolver.fvSchemes.divSchemes['div\(phi.*,U.*\)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,U)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = 'Gauss linear'
thSolver.fvSchemes.divSchemes['div(alphaRhoPhi,K)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,k.*\)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,epsilon.*\)'] = 'Gauss upwind'
thSolver.fvSchemes.divSchemes['div\(alphaRhoPhi.*,(h|e).*\)'] = 'Gauss upwind'

thSolver.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

thSolver.fvSchemes.snGradSchemes['default'] = 'uncorrected'


#==============================================================================*
# Neutronics solver

neutronicsSolver = ffn.solvers.NeutronicsSolver(
    region=nMesh.region,
    solver="diffusionNeutronics",
    mesh=nMesh,
    power=10000000,
    fastNeutrons=True,
)

pointKineticsData = ffn.nuclearData.PointKineticsData(
    promptGenerationTime=1e-06,
    delayedFractions=[
        7.2315e-05,
        0.000609661,
        0.000471181,
        0.00118907,
        0.000445487,
        9.58515e-05,
    ],
    decayConstants=[
        0.0125371,
        0.0300828,
        0.109879,
        0.325484,
        1.3036,
        9.51817,
    ],
    feedbackCoeffTFuel=-3e-06
)

pointKineticsData.controlRodReactivityMap = [
    (  0.1, -0.01 ),
    (  0.0,  0.0  ),
    ( -0.1,  0.01 ),
]

pointKineticsData.externalReactivityTimeProfile = ffn.TimeProfile(
    'table',
    startTime=0,
    table=[
        (   0,   0.0005767131    )
    ]
)

# neutronicsSolver.nuclearData = pointKineticsData
neutronicsSolver.nuclearData.add_state(
    ffn.nuclearData.NuclearDataState(
        name='reference',
        zones=[
            ffn.nuclearData.NuclearDataZone(
                name='zone0',
                fuelFraction=1,
                inverseVelocity=[7.1e-08, 0.0045],
                diffusionCoefficient=[0.1, 0.04],
                nuFissionXS=[1, 1.75],
                powerXS=[1, 1],
                scatteringMatrixP0=[[10, 5], [0, 10]],
                removalXS=[5, 4],
                chiPrompt=[1, 0],
                chiDelayed=[0.5, 0.5],
                delayedFraction=[7.2315e-05, 0.000609661, 0.000471181, 0.00118907, 0.000445487, 9.58515e-05],
                decayConstant=[0.0125371, 0.0300828, 0.109879, 0.325484, 1.3036, 9.51817],
                discFactor=[1, 1],
                integralFlux=[1, 1]
            )
        ]
    )
)

neutronicsSolver.fvSchemes.divSchemes['div(facePhi_,angularFlux_)'] = "Gauss upwind"

neutronicsSolver.neutronTransportOptions.integralPredictor = True
neutronicsSolver.neutronTransportOptions.aitkenAcceleration = True
neutronicsSolver.neutronTransportOptions.maxNeutronIterations = 50

#==============================================================================*
# Solvers

solvers = ffn.solvers.Solvers([neutronicsSolver, thSolver])

#==============================================================================*
# Coupling

coupling = ffn.coupling.Coupling(solvers)

coupling.add_field_transfer(neutronicsSolver, thSolver, 'powerDensity', 'powerDensityStructure')
coupling.add_field_transfer(neutronicsSolver, thSolver, 'secondaryPowerDensity', 'powerDensityLiquid')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T', 'TCool')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'thermo:rho', 'rhoCool')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.fuelAvForNeutronics', 'TFuel')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.cladAvForNeutronics', 'TClad')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.passiveStructure', 'TStructMech')


#==============================================================================*
# Model

model = ffn.Case(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0]
)

settings = model.settings

settings.application = "GeN-Foam"
settings.startTime = 0
settings.endTime = 100
settings.deltaT = 1
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = 100
settings.adjustTimeStep = True

print(model)

model.coupling.plot_coupling_graph()
model.coupling.plot_solving_graph()
model.coupling.plot_solving_flowchart()

model.export_to_openfoam()


#==============================================================================*
# Run

ffn.run(model, is_preprocessing=True)


#==============================================================================*
# Post-processing

model.plot_slice(
    region=nMesh.region,
    time=settings.endTime,
    fieldName='powerDensity',
    offset=(0, 0.05, 0)
)
model.plot_slice(
    region=thMesh.region,
    time=settings.endTime,
    fieldName='T',
    offset=(0, 0.05, 0),
    cmap='RdBu_r',
    show_edges=True,
    unit='K'
)

model.plot_residuals(
    parameters=['fluxStar0'],
    title="Neutronics"
)
model.plot_residuals(
    parameters=['p_rgh', 'h'],
    title="Thermal-hydraulics"
)

#==============================================================================*
