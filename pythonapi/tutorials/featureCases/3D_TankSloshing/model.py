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

(centerBlock, frontBlock, rightBlock, backBlock, leftBlock) = thMesh.create_cylinder_along_z(
    'tank', 1, -tankLength/2, tankLength/2, 10, 10, 30
)

wall = mesh.Face("wall", boundaryType="wall")
for block in [centerBlock, frontBlock, rightBlock, backBlock, leftBlock]:
    wall.add_sub_face(block.topFace())
    wall.add_sub_face(block.bottomFace())

wall.add_sub_face(frontBlock.frontFace())
wall.add_sub_face(rightBlock.rightFace())
wall.add_sub_face(backBlock.backFace())
wall.add_sub_face(leftBlock.leftFace())

thMesh.add_boundary(wall)


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

U = ffn.fields.Field('U', region=thMesh.region)
U.dimensions = ffn.fields.Dimension('U')
U.internalField = ffn.common.Vector(0, 0, 0)
U.set_boundary_condition("wall", bc.MovingWallVelocity(U.internalField))

T = ffn.fields.Field('T', region=thMesh.region)
T.dimensions = ffn.fields.Dimension('T')
T.internalField = 300
T.set_boundary_condition("wall", bc.FixedValue(T.internalField))

p = ffn.fields.Field("p", region=thMesh.region)
p.dimensions = ffn.fields.Dimension(default='p')
p.internalField = 1.1e5
p.set_boundary_condition("wall", bc.FixedFluxPressure(p.internalField))

p_rgh = ffn.fields.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.fields.Dimension(default='p')
p_rgh.internalField = 1.1e5
p_rgh.set_boundary_condition("wall", bc.FixedFluxPressure(p_rgh.internalField))

alpha_water = ffn.fields.Field("alpha.water", region=thMesh.region)
alpha_water.dimensions = ffn.fields.Dimension()
alpha_water.internalField = 0
alpha_water.set_boundary_condition("wall", bc.ZeroGradient())

timeFolder0.append(U)
timeFolder0.append(T)
timeFolder0.append(p)
timeFolder0.append(p_rgh)
timeFolder0.append(alpha_water)

#==============================================================================*
# Solver

thSolver = ffn.solvers.thermal_hydraulics.CompressibleInterFoam(
    region=thMesh.region,
    mesh=thMesh
)

thSolver.fluid.turbulenceProperties.simulationType = 'laminar'

multiPhase = ffn.thermo.MultiPhase()

multiPhase.appendPhase(ffn.thermo.WaterPerfectFluid(ext='water'))
multiPhase.appendPhase(ffn.thermo.AirPerfectGas(ext='air'))

thSolver.fluid.thermophysicalProperties = multiPhase

thSolver.setFieldsDict.add_default_values(alpha_water, 0)


# Filling the tank
#-----------------

thSolver.setFieldsDict.add_box_to_cell(
    fieldValues=[(alpha_water, 1)],
    lowCorner=ffn.common.Vector(-10, -10, -10),
    highCorner=ffn.common.Vector(10, 10, 0)
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


# Motion setup
#-------------

thSolver.dynamicMeshDict.dynamicFvMesh = "dynamicMotionSolverFvMesh"
thSolver.dynamicMeshDict.motionSolver = "solidBody"
thSolver.dynamicMeshDict.solidBodyMotionFunction = "multiMotion"
omegaTable = ffn.common.Table([
    (3, 0.5*np.pi),
    (7, -0.5*np.pi),
])

# Choose one or more of the following option

# --- Linear motion
# thSolver.dynamicMeshDict.add_linear_motion(name='motion1', velocity=ffn.common.Vector(0.2, 0, 0))

# --- Linear oscillation long an axis
# thSolver.dynamicMeshDict.add_oscillating_linear_motion('motion21', ffn.common.Vector(0.5, 0, 0), omega=0.5*np.pi)
# thSolver.dynamicMeshDict.add_oscillating_linear_motion('motion21', ffn.common.Vector(0.5, 0, 0), omega=omegaTable)
# thSolver.dynamicMeshDict.add_oscillating_linear_motion('motion22', ffn.common.Vector(0, 1.0, 0), omega=1.0*np.pi)

# --- Rotation motion
# thSolver.dynamicMeshDict.add_rotating_motion('motion3', 0.5*np.pi, axis=ffn.common.Vector(1, 0, 0))

# --- Rotation oscillation
# thSolver.dynamicMeshDict.add_oscillating_rotating_motion('motion4', 0.5*np.pi, amplitude=ffn.common.Vector(0, 45, 0))

# --- User-defined position and rotation over time
dofMotion = mesh.dicts.Tabulated6DoFMotion(
    'tank',
    tableDoF=[
        (0,   ffn.common.Vector(0, 0, 0)),
        (2.5, ffn.common.Vector(0, -0.5, -0.5)),
        (5,   ffn.common.Vector(0, 0.5, 0), ffn.common.Vector(0, 90, 0)),
        (7.5, ffn.common.Vector(0, 0.5, 0.5)),
        (10,  ffn.common.Vector(0, 0, 0)),
    ]
)
# or
# dofMotion.add_point(0,   ffn.common.Vector(0, 0, 0))
# dofMotion.add_point(2.5, ffn.common.Vector(0, -0.5, -0.5))
# dofMotion.add_point(5,   ffn.common.Vector(0, 0.5, 0), ffn.common.Vector(0, 90, 0))
# dofMotion.add_point(7.5, ffn.common.Vector(0, 0.5, 0.5))
# dofMotion.add_point(10,  ffn.common.Vector(0, 0, 0))

thSolver.dynamicMeshDict.add_motion(dofMotion)


#==============================================================================*
# Model

solvers = ffn.solvers.Solvers([thSolver])

model = ffn.Case(
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
