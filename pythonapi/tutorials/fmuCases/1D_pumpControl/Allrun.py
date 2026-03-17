"""

"""
#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh
import numpy as np


#==============================================================================*
# Mesh

pipeWidth = 0.5

thMesh = mesh.BlockMesh(region="fluidRegion")

pipe1 = thMesh.create_cube("pipe1", 0, -pipeWidth/2, -pipeWidth/2, 1, pipeWidth/2, pipeWidth/2, nx=10)
pump = thMesh.extrude_right([pipe1], "pump", dx=0.2, nx=2)
pipe2 = thMesh.extrude_right([pump], "pipe2", dx=1, nx=10)

pipe3 = thMesh.add_pipe_1D_from_direction(
    "pipe3",
    originPosition=pipe2,
    direction=ffn.common.Vector(0, 0, 1),
    length=2,
    n=10,
    equivalentHydraulicDiameter=pipeWidth/np.sqrt(np.pi),
    elbowRadius=pipeWidth,
    isAddBoundaryConditions=True,
    originPositionOutletFaceName='right'
)
pipe4 = thMesh.add_pipe_1D_from_direction(
    "pipe4",
    originPosition=pipe3,
    direction=ffn.common.Vector(-1, 0, 0),
    length=2.2,
    n=10,
    equivalentHydraulicDiameter=pipeWidth/np.sqrt(np.pi),
    elbowRadius=pipeWidth,
    isAddBoundaryConditions=True,
)
pipe5 = thMesh.add_pipe_1D_from_2points(
    "pipe5",
    originPosition=pipe4,
    finalPosition=pipe1,
    n=10,
    equivalentHydraulicDiameter=pipeWidth/np.sqrt(np.pi),
    elbowRadius=pipeWidth,
    isAddBoundaryConditions=True,
    finalPositionInletFaceName='left'
)

walls = mesh.Face("walls", boundaryType="wall")
for block in [pipe1, pump, pipe2]:
    walls.add_sub_face(block.frontFace())
    walls.add_sub_face(block.backFace())
    walls.add_sub_face(block.topFace())
    walls.add_sub_face(block.bottomFace())

thMesh.add_boundary(walls)


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

T = ffn.fields.Field("T", region=thMesh.region)
T.dimensions = ffn.fields.Dimension(default='T')
T.internalField = 600
T.set_boundary_condition(walls, bc.ZeroGradient())

U = ffn.fields.Field("U", region=thMesh.region)
U.dimensions = ffn.fields.Dimension(default='U')
U.internalField = ffn.common.Vector(0, 0, 0)
U.set_boundary_condition(walls, bc.Slip())

p_rgh = ffn.fields.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.fields.Dimension(default='p')
p_rgh.internalField = 155e5
p_rgh.set_boundary_condition(walls, bc.ZeroGradient())

for field in [p_rgh, U, T]:
    field.set_boundary_condition('"pipe.*_outlet"', bc.Cyclic())
    field.set_boundary_condition('"pipe.*_inlet"', bc.Cyclic())
for field in [p_rgh, T]:
    field.set_boundary_condition(thMesh.pipeWallBC, bc.ZeroGradient())
U.set_boundary_condition(thMesh.pipeWallBC, bc.Slip())

timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p_rgh)


#==============================================================================*
# Thermal-hydraulics solver

thSolver = ffn.solvers.thermal_hydraulics.OnePhaseThermalHydraulicsSolver(
    region=thMesh.region,
    mesh=thMesh,
    removeBaffles=True,
    isSetFvSolutionToDefault=False
)

thSolver.fluid.thermophysicalProperties = ffn.thermo.WaterPolynomial()

thSolver.fluid.turbulenceProperties.simulationType = 'laminar'

pumpModel = ffn.porous_medium.Pump(
    zones=['pump'],
    volumeFraction=0.2,
    Dh=0.1,
    momentumSource=ffn.common.Vector(1e7, 0, 0),
    momentumSourceTimeProfile=ffn.TimeProfile(
        type='fmi',
        nameFromFMU="gfMomentumSource",
        initialValue=0
    )
)

flowBlockageModel = ffn.porous_medium.Structure(
    zones=['pipe3'],
    pitch=0.02,
    elementDiameter=0.007,
    latticeType='hexagon'
)

thSolver.structures.append(pumpModel)
thSolver.structures.append(flowBlockageModel)


thSolution = ffn.numerics.fvSolution()

thSolution.append('"p_rgh.*"', ffn.numerics.fvSolutionSolver(
    solver='GAMG', smoother='DIC', tolerance=1e-8, relTol=0
))
thSolution.append('"e.*"', ffn.numerics.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.append('"h.*"', ffn.numerics.fvSolutionSolver(
    solver='smoothSolver',
    smoother='symGaussSeidel',
    tolerance=1e-8, relTol=0, minIter=0
))
thSolution.append('".*"', ffn.numerics.fvSolutionSolver(
    solver='PBiCGStab',
    preconditioner='diagonal',
    tolerance=1e-6, relTol=0.001
))

thSolver.fvSolution = thSolution

thSolver.pimpleOptions.nOuterCorrectors = 3
thSolver.pimpleOptions.nCorrectors = 3
thSolver.pimpleOptions.pMin = 150e5


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

thSolver.add_relaxation_on_field('".*"', 0.7)
thSolver.add_relaxation_on_field('"p_rgh.*"', 1)
thSolver.add_relaxation_on_equation('".*"', 1)
thSolver.add_relaxation_on_equation('"e.*"', 0.7)
thSolver.add_relaxation_on_equation('"k.*"', 0.7)
thSolver.add_relaxation_on_equation('"h.*"', 0.7)


#==============================================================================*
# Solvers

solvers = ffn.solvers.Solvers([thSolver])


#==============================================================================*
# FMU
# The FMU needs to be linked even during steady-state, this improves stability
# of the results especially if PID controllers are present

# Coupling interface
externalCouplingDict = ffn.ExternalCouplingDict()

massFlowRateSens = ffn.MassFlowToFMU(
    name="massFlowRateSens",
    nameFMU='gfMassFlow',
    regionName='pipe2_outlet',
    regionType='patch',
    region=thMesh.region,
    initValue=300
)
externalCouplingDict.append(massFlowRateSens)

# inletPsensor = ffn.ExtSensor(
#     name='inletP',
#     fieldName='p_rgh',
#     sensorName="gfInletP",
#     sensorPosition=ffn.common.Vector(0.9, pipeWidth/2, pipeWidth/2),
#     region=thMesh.region,
#     initValue=155e5
# )
# outletPsensor = ffn.ExtSensor(
#     name='outletP',
#     fieldName='p_rgh',
#     sensorName="gfOutletP",
#     sensorPosition=ffn.common.Vector(1.3, pipeWidth/2, pipeWidth/2),
#     region=thMesh.region,
#     initValue=155e5
# )
# externalCouplingDict.append(inletPsensor)
# externalCouplingDict.append(outletPsensor)

# FMU simulator
FMUSimulator = ffn.FMUSimulator(
    name="FMUSimulator",
    pyClassName="Controller",
    pyFileName="Controller",
    mapping=[
        ("from", "gfMassFlow", "massFlow"),
        # ("from", "gfInletP", "inletP"),
        # ("from", "gfOutletP", "outletP"),
        ("to", "gfMomentumSource", "pumpMomentum")
    ]
)
FMUSimulator.plot_coupling_graph()


#==============================================================================*
# Model

model = ffn.Case(
    solvers=solvers,
    timeFolders=[timeFolder0],
    externalCouplingDict=externalCouplingDict
)

model.add_function_object(FMUSimulator)

settings = model.settings

settings.application = "GeN-Foam"
settings.startTime = 0
settings.endTime = 200
settings.deltaT = 0.1
settings.startFrom = 'latestTime'
settings.writeControl = 'adjustableRunTime'
settings.writeInterval = settings.endTime/10
settings.adjustTimeStep = True
settings.maxDeltaT = 0.1
settings.maxCo = 1
settings.maxPowerVariation = 0.01

print(model)

ffn.allclean()
model.export_to_openfoam()


#==============================================================================*
# Run

ffn.run_preprocessing(model=model)

model.plot_mesh(region=thMesh, show_edges=True)
model.plot_slice(region=thMesh, show_edges=True)

ffn.run(model, is_preprocessing=False)


#==============================================================================*
# Post-processing

model.plot_mesh(
    region=thMesh,
    time=model.settings.endTime,
    fieldName='U',
    cmap='RdBu_r',
    normal='y',
    show_edges=True,
    unit='m/s'
)
model.plot_mesh(
    region=thMesh,
    time=model.settings.endTime,
    fieldName='p_rgh',
    cmap='RdBu_r',
    normal='y',
    show_edges=True,
    unit='Pa'
)

import Allpostprocess

#==============================================================================*