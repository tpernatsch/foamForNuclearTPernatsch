"""

"""
#==============================================================================*
# Imports

from matplotlib import pyplot as plt
import numpy as np
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# Mesh

tankLength = 3

thMesh = mesh.BlockMesh(region="")

(centerBlock, frontBlock, rightBlock, backBlock, leftBlock) = thMesh.createCylinderAlongZ(
    'tank', 1, -tankLength/2, tankLength/2, 10, 10, 30
)

wall = ffn.Face("wall", boundaryType="wall")
for block in [centerBlock, frontBlock, rightBlock, backBlock, leftBlock]:
    wall.addSubFace(block.topFace())
    wall.addSubFace(block.bottomFace())

wall.addSubFace(frontBlock.frontFace())
wall.addSubFace(rightBlock.rightFace())
wall.addSubFace(backBlock.backFace())
wall.addSubFace(leftBlock.leftFace())

thMesh.addBoundary(wall)


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

U = ffn.Field('U', region=thMesh.region)
U.dimensions = ffn.Dimension('U')
U.internalField = ffn.Vector(0, 0, 0)
U.set_boundary_condition("wall", bc.MovingWallVelocity(U.internalField))

T = ffn.Field('T', region=thMesh.region)
T.dimensions = ffn.Dimension('T')
T.internalField = 300
T.set_boundary_condition("wall", bc.FixedValue(T.internalField))

p = ffn.Field("p", region=thMesh.region)
p.dimensions = ffn.Dimension(default='p')
p.internalField = 1.1e5
p.set_boundary_condition("wall", bc.FixedFluxPressure(p.internalField))

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 1.1e5
p_rgh.set_boundary_condition("wall", bc.FixedFluxPressure(p_rgh.internalField))

alpha_water = ffn.Field("alpha.water", region=thMesh.region)
alpha_water.dimensions = ffn.Dimension()
alpha_water.internalField = 0
alpha_water.set_boundary_condition("wall", bc.ZeroGradient())

timeFolder0.append(U)
timeFolder0.append(T)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(alpha_water)

#==============================================================================*
# Solver

thSolver = ffn.CompressibleInterFoam(region=thMesh.region, mesh=thMesh)

thSolver.turbulenceProperties.simulationType = 'laminar'

multiPhase = ffn.thermophysicalProperty.MultiPhase()

multiPhase.appendPhase(ffn.thermophysicalProperty.WaterPerfectFluid(ext='water'))
multiPhase.appendPhase(ffn.thermophysicalProperty.AirPerfectGas(ext='air'))

thSolver.thermophysicalProperties = multiPhase

thSolver.setFieldsDict.add_default_values(alpha_water, 0)

# Fill half the tank
thSolver.setFieldsDict.add_box_to_cell(
    fieldValues=[(alpha_water, 1)],
    lowCorner=ffn.Vector(-10, -10, -10),
    highCorner=ffn.Vector(10, 10, 0)
)
# thSolver.setFieldsDict.add_cylinder_annulus_to_cell(
#     fieldValues=[(alpha_water, 1)],
#     innerRadius=0.5,
#     outerRadius=0.8,
#     point1=ffn.Vector(0, 0, -0.4),
#     point2=ffn.Vector(0, 0, 0.4),
# )
# thSolver.setFieldsDict.add_sphere_to_cell(
#     fieldValues=[(alpha_water, 1)],
#     outerRadius=0.4,
#     innerRadius=0.1,
#     origin=ffn.Vector(0, 0, 0.5)
# )

thSolver.dynamicMeshDict.dynamicFvMesh = "dynamicMotionSolverFvMesh"
thSolver.dynamicMeshDict.motionSolver = "solidBody"
thSolver.dynamicMeshDict.solidBodyMotionFunction = "multiMotion"
omegaTable = ffn.Table([
    (3, 0.5*np.pi),
    (7, -0.5*np.pi),
])
# thSolver.dynamicMeshDict.add_linear_motion(name='motion1', velocity=ffn.Vector(0.2, 0, 0))
thSolver.dynamicMeshDict.add_oscillating_linear_motion('motion21', ffn.Vector(0.5, 0, 0), omega=0.5*np.pi)
# thSolver.dynamicMeshDict.add_oscillating_linear_motion('motion21', ffn.Vector(0.5, 0, 0), omega=omegaTable)
# thSolver.dynamicMeshDict.add_oscillating_linear_motion('motion22', ffn.Vector(0, 1.0, 0), omega=1.0*np.pi)
# thSolver.dynamicMeshDict.add_rotating_motion('motion3', 0.5*np.pi, axis=ffn.Vector(1, 0, 0))
# thSolver.dynamicMeshDict.add_oscillating_rotating_motion('motion4', 0.5*np.pi, amplitude=ffn.Vector(0, 45, 0))


# thSolver.dynamicMeshDict.dynamicFvMesh = "dynamicMotionSolverFvMesh"
# thSolver.dynamicMeshDict.motionSolver = "velocityLaplacian"
# thSolver.dynamicMeshDict.diffusivity = ffn.QuadraticDiffusivityMotion('inverseDistance', ['wall'])

#==============================================================================*
# Model

solvers = ffn.Solvers([thSolver])

model = ffn.Model(
    solvers=solvers,
    timeFolders=[timeFolder0]
)
model.settings.application = "compressibleInterDyMFoam"
model.settings.endTime = 10
model.settings.deltaT = 0.001
model.settings.writeControl = "adjustableRunTime"
model.settings.writeInterval = 0.1
model.settings.adjustTimeStep = True
model.settings.runTimeModifiable = True
model.settings.maxCo = 1
model.settings.maxDeltaT = 1
model.settings.maxAlphaCo = 1

fig, ax = plt.subplots(figsize=(5, 4), dpi=200)
tRange = np.linspace(model.settings.startTime, model.settings.endTime, 100)
ax.plot(tRange, [omegaTable.value(t) for t in tRange])
ax.set_xlabel("Time [s]")
ax.set_ylabel("Omega [rad/s]")
fig.tight_layout()
fig.savefig("fig_data_omegaTable.png")

#==============================================================================*
