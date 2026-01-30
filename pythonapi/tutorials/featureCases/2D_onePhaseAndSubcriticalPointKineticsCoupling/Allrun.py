"""

"""
#==============================================================================*
# Imports

from copy import copy
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*

ffn.allclean()

#==============================================================================*
# Mesh

nMesh = mesh.BlockMesh(region="neutroRegion")

zone0 = nMesh.create_cube("zone0", 0, 0, 0, 1, 0.1, 1, nx=25, nz=25)

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

nMesh.add_boundary(inlet)
nMesh.add_boundary(outlet)
nMesh.add_boundary(fixedWalls)
nMesh.add_boundary(frontAndBack)


thMesh = copy(nMesh)
thMesh.region = "fluidRegion"


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

defaultFlux = ffn.Field("defaultFlux", region=nMesh.region)
defaultFlux.dimensions = ffn.Dimension(default='neutronFlux')
defaultFlux.internalField = 1
defaultFlux.set_boundary_condition(fixedWalls, bc.FixedValue(0))
defaultFlux.set_boundary_condition(inlet, bc.FixedValue(0))
defaultFlux.set_boundary_condition(outlet, bc.FixedValue(0))
defaultFlux.set_boundary_condition(frontAndBack, bc.Empty())

T = ffn.Field("T", region=thMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = 600
T.set_boundary_condition(fixedWalls, bc.ZeroGradient())
T.set_boundary_condition(inlet, bc.FixedValue(600))
T.set_boundary_condition(outlet, bc.ZeroGradient())
T.set_boundary_condition(frontAndBack, bc.Empty())

U = ffn.Field("U", region=thMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 1)
U.set_boundary_condition(fixedWalls, bc.Slip())
U.set_boundary_condition(inlet, bc.FixedValue(U.internalField))
U.set_boundary_condition(outlet, bc.ZeroGradient())
U.set_boundary_condition(frontAndBack, bc.Empty())

p = ffn.Field("p", region=thMesh.region)
p.dimensions = ffn.Dimension(default='p')
p.internalField = 1e5
p.set_boundary_condition(fixedWalls, bc.Calculated(p.internalField))
p.set_boundary_condition(inlet, bc.Calculated(p.internalField))
p.set_boundary_condition(outlet, bc.Calculated(p.internalField))
p.set_boundary_condition(frontAndBack, bc.Empty())

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 1e5
p_rgh.set_boundary_condition(fixedWalls, bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition(inlet, bc.FixedFluxPressure(p_rgh.internalField))
p_rgh.set_boundary_condition(outlet, bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition(frontAndBack, bc.Empty())

timeFolder0.append(defaultFlux)
timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p)
timeFolder0.append(p_rgh)


#==============================================================================*
# Thermal-hydraulics solver

thSolver = ffn.ThermalHydraulicsSolver(
    "fluidRegion", "onePhase",
    mesh=thMesh,
    removeBaffles=True,
    isSetFvSolutionToDefault=False
)

thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumConst()

thSolver.turbulenceProperties.simulationType = 'laminar'

core = ffn.StructureProperty(['zone0'], volumeFraction=0.5, Dh=0.01)

core.powerModel = ffn.NuclearFuelPin(
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

core.add_passive_structure(
    volumetricArea=2,
    rho=7700,
    Cp=500,
    T=600
)

thSolver.add_structure_property(core)

thSolver.add_drag_model(
    ffn.ReynoldsPower(0.687, -0.25, zones=['zone0'])
)
thSolver.add_heat_transfer_model(
    ffn.NusseltReynoldsPrandtlPower(4.82, 0.0185, 0.827, 0.827, zones=['zone0'])
)


thSolution = ffn.fvSolution()

thSolution.add_fvSolutionSolver('p_rgh.*', ffn.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-8, relTol=0
))
thSolution.add_fvSolutionSolver('e.*', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.add_fvSolutionSolver('h.*', ffn.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.add_fvSolutionSolver('.*', ffn.fvSolutionSolver(
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

neutronicsSolver = ffn.NeutronicsSolver(
    "neutroRegion",
    "diffusionNeutronics",
    mesh=nMesh,
    power=10000000,
    keff=1,
    ksrc=0.98
)

neutronicsSolver.nuclearData.add_state(
    ffn.NuclearDataState(
        name='reference',
        zones=[
            ffn.NuclearDataZone(
                name='zone0',
                fuelFraction=1,
                inverseVelocity=[7.1e-08, 0.0045],
                diffusionCoefficient=[0.1, 0.04],
                nuFissionXS=[1, 5.74],
                powerXS=[1, 1],
                scatteringMatrixP0=[[10, 5], [0, 10]],
                removalXS=[5, 4],
                chiPrompt=[1, 0],
                chiDelayed=[0.5, 0.5],
                delayedFraction=[0.0075],
                decayConstant=[1e-5],
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

solvers = ffn.Solvers([neutronicsSolver, thSolver])

#==============================================================================*
# Coupling

coupling = ffn.Coupling(solvers)

coupling.add_field_transfer(neutronicsSolver, thSolver, 'powerDensity', 'powerDensityStructure')
coupling.add_field_transfer(neutronicsSolver, thSolver, 'secondaryPowerDensity', 'powerDensityLiquid')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T', 'TCool')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'thermo:rho', 'rhoCool')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.fuelAvForNeutronics', 'TFuel')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.cladAvForNeutronics', 'TClad')
coupling.add_field_transfer(thSolver, neutronicsSolver, 'T.passiveStructure', 'TStructMech')


#==============================================================================*
# Model

transientStartTime = 150

model = ffn.Model(
    solvers=solvers,
    coupling=coupling,
    timeFolders=[timeFolder0],
    caseFolder="steadyState"
)

settings = model.settings

settings.application = "GeN-Foam"
settings.startTime = 0
settings.endTime = transientStartTime
settings.deltaT = 1
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = transientStartTime
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
    fieldName='powerDensity'
)
model.plot_slice(
    region=thMesh.region,
    time=settings.endTime,
    fieldName='T',
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
# Transient Source

neutronicsSolver.solver = "pointKinetics"
neutronicsSolver.eigenvalueNeutronics = False
neutronicsSolver.externalSourceNeutronics = True

neutronicsSolver.externalSource.nuSource = 15.5 # [neutron/proton]
neutronicsSolver.externalSource.beamEnergy = 1.28160e-10 # [J/proton], 800 MeV/proton

neutronicsSolver.externalSource.externalSourceModulationTimeProfile = ffn.TimeProfile(
    type="table",
    startTime=transientStartTime,
    table=[
        (0,  1),
        (20, 2)
    ]
)

pointKineticsData = ffn.PointKineticsData(
    fastNeutrons=True,
    promptGenerationTime=0.0015,
    nuFission=2.45,
    energyPerFission=32e-12,
    delayedFractions=[0.000075],
    decayConstants=[1e-5],
    # feedbackCoeffTFuel=-3e-06
)
pointKineticsData.externalReactivityTimeProfile = ffn.TimeProfile(
    'table',
    startTime=1000000000,
    table=[
        (0,    0),
        (0.01, 1020.408163265307e-5)
    ]
)

neutronicsSolver.nuclearData = pointKineticsData

model.caseFolder = "transientSource"

settings.endTime = settings.endTime + 100
settings.deltaT = 1e-7
settings.writeControl = "adjustableRunTime"
settings.writeInterval = 10

ffn.copyFolder("steadyState", model.caseFolder)

print(model)

model.export_to_openfoam()

ffn.run(model=model, is_preprocessing=False)


#==============================================================================*
# Transient Reactivity

neutronicsSolver.externalSource.externalSourceModulationTimeProfile.startTime = 1e32
pointKineticsData.externalReactivityTimeProfile.startTime = transientStartTime

model.caseFolder = "transientReactivity"

ffn.copyFolder("steadyState", model.caseFolder)

model.export_to_openfoam()

ffn.run(model=model, is_preprocessing=False)


#==============================================================================*
