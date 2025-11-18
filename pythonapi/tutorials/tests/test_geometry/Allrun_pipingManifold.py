#==============================================================================*
# Imports

import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh
import foamForNuclear.boundaryConditions as bc


#==============================================================================*
# Mesh

thMesh = mesh.BlockMesh(region='fluidRegion')

coreRadius = 2
coreLength = 3
legLength = 3

theta = 30 * np.pi/180

elbowRadius = 0.5
equivalentHydraulicDiameter = 0.5

core = thMesh.create_cylinder_along_z(
    name='core',
    radius=coreRadius,
    lowZ=0, highZ=coreLength,
    nx=3, ny=3, nz=10,
    isAddAllBC=True
)

topSections = thMesh.create_pipe_cylindrical_manifold_along_z(
    'manifold',
    nEntries=2,
    innerRadius=coreRadius/2,
    outerRadius=coreRadius,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    lowZ=coreLength,
    nr=3, nt=6
)

manifoldInlet = ffn.Face("manifoldInlet")
for section in topSections:
    manifoldInlet.add_sub_face(section.bottomFace())
thMesh.add_boundary(manifoldInlet)

thMesh.merge_patch_pairs_by_name("coreTop_1", manifoldInlet.name)


manifoldPipes = topSections[::2]

hotLeg1 = thMesh.extrude_normal(manifoldPipes[0], 'front', 'hotLeg1', 1, 4)

hotLeg3 = thMesh.add_pipe_1D_from_direction(
    "hotLeg3",
    originPosition=ffn.Vector(4, -1, coreLength+np.sqrt(np.pi)*equivalentHydraulicDiameter/2),
    direction=ffn.Vector(0, -1, 0),
    length=legLength,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    elbowRadius=elbowRadius,
    n=10,
    isAddLateralBC=True
)

hotLeg2 = thMesh.add_pipe_1D_from_2points(
    "hotLeg2",
    originPosition=hotLeg1,
    finalPosition=hotLeg3,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    elbowRadius=elbowRadius,
    n=4,
    isAddLateralBC=True,
    originPositionOutletFaceName='front'
)

thMesh.merge_patches_with_name(name="walls", includeFacename=["coreWall", "coreTop"], patchType="wall")
thMesh.merge_patches_with_name(name="coreInlet", includeFacename=["coreBottom"])

hotLegOutlet = ffn.Face("hotLegOutlet")
hotLegOutlet.add_sub_face(hotLeg3.topFace())
thMesh.add_boundary(hotLegOutlet)

# thMesh.emptyBoundaryCondition.boundaryType = 'patch'


#==============================================================================*
# Time folder

timeFolder0 = ffn.TimeFolder(0)

T = ffn.Field("T", region=thMesh.region)
T.dimensions = ffn.Dimension(default='T')
T.internalField = 600
T.set_boundary_condition("pipeWall", bc.ZeroGradient())
T.set_boundary_condition("hotLegOutlet", bc.ZeroGradient())
T.set_boundary_condition("coreInlet", bc.FixedValue(T.internalField))
T.set_boundary_condition("walls", bc.ZeroGradient())
T.set_boundary_condition("manifoldInlet", bc.ZeroGradient())
T.set_boundary_condition("coreTop_1", bc.ZeroGradient())
T.set_boundary_condition("defaultFaces", bc.Empty())

U = ffn.Field("U", region=thMesh.region)
U.dimensions = ffn.Dimension(default='U')
U.internalField = ffn.Vector(0, 0, 0)
U.set_boundary_condition("pipeWall", bc.Slip())
U.set_boundary_condition("hotLegOutlet", bc.ZeroGradient())
U.set_boundary_condition("coreInlet", bc.ZeroGradient())
U.set_boundary_condition("walls", bc.Slip())
U.set_boundary_condition("manifoldInlet", bc.ZeroGradient())
U.set_boundary_condition("coreTop_1", bc.Slip())
U.set_boundary_condition("defaultFaces", bc.Empty())

p_rgh = ffn.Field("p_rgh", region=thMesh.region)
p_rgh.dimensions = ffn.Dimension(default='p')
p_rgh.internalField = 1e5
p_rgh.set_boundary_condition("pipeWall", bc.ZeroGradient())
p_rgh.set_boundary_condition("hotLegOutlet", bc.FixedValue(p_rgh.internalField))
p_rgh.set_boundary_condition("coreInlet", bc.FixedValue(3*p_rgh.internalField))
p_rgh.set_boundary_condition("walls", bc.ZeroGradient())
p_rgh.set_boundary_condition("manifoldInlet", bc.ZeroGradient())
p_rgh.set_boundary_condition("coreTop_1", bc.ZeroGradient())
p_rgh.set_boundary_condition("defaultFaces", bc.Empty())

for face in [hotLeg1, hotLeg2]:
    T.set_boundary_condition(face.name+"_outlet", bc.Cyclic())
    U.set_boundary_condition(face.name+"_outlet", bc.Cyclic())
    p_rgh.set_boundary_condition(face.name+"_outlet", bc.Cyclic())
    # p_rgh.set_boundary_condition(face, bc.FixedJump(value=p_rgh.internalField, jump=0, rho='"thermo:rho"'))

for face in [hotLeg2, hotLeg3]:
    T.set_boundary_condition(face.name+"_inlet", bc.Cyclic())
    U.set_boundary_condition(face.name+"_inlet", bc.Cyclic())
    p_rgh.set_boundary_condition(face.name+"_inlet", bc.Cyclic())
    # p_rgh.set_boundary_condition(face, bc.FixedJump(value=p_rgh.internalField, jump=0, rho='"thermo:rho"'))

timeFolder0.append(T)
timeFolder0.append(U)
timeFolder0.append(p_rgh)


#==============================================================================*
# Solver

thSolver = ffn.ThermalHydraulicsSolver(
    region=thMesh.region,
    mesh=thMesh,
    solver='onePhase',
    isSetFvSchemesToDefault=False
)

thSolver.thermophysicalProperties = ffn.thermophysicalProperty.SodiumConst()

thSolver.turbulenceProperties.simulationType = 'laminar'

corePowerModel = ffn.StructureProperty(['core'], volumeFraction=0.5, Dh=0.01)

corePowerModel.powerModel = ffn.NuclearFuelPin(
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
    fuelT=T.internalField,
    cladT=T.internalField,
    gapH=3000,
    powerDensity=10e6
)

corePowerModel.add_passive_structure(
    volumetricArea=2,
    rho=7700,
    Cp=500,
    T=T.internalField
)

thSolver.add_structure_property(corePowerModel)

thSolver.add_drag_model(
    ffn.ReynoldsPower(0.687, -0.25, zones=['core'])
)
thSolver.add_heat_transfer_model(
    ffn.NusseltReynoldsPrandtlPower(4.82, 0.0185, 0.827, 0.827, zones=['core'])
)

thSolver.pimpleOptions.nOuterCorrectors = 2
thSolver.pimpleOptions.nCorrectors = 1
thSolver.pimpleOptions.nNonOrthogonalCorrectors = 0
thSolver.pimpleOptions.correctPhi = False
thSolver.pimpleOptions.solveEnergy = True

thSolver.fvSchemes.divSchemes['default'] = "Gauss upwind"
thSolver.fvSchemes.divSchemes['div(alphaRhoPhiNu,U)'] = "Gauss linear"


#==============================================================================*
# Model

model = ffn.Model(timeFolders=[timeFolder0])
model.settings.application = "GeN-Foam"
model.settings.endTime = 1000
model.settings.deltaT = 0.001
model.settings.writeInterval = 10
model.settings.writeControl = "adjustableRunTime"
model.settings.adjustTimeStep = True
model.settings.runTimeModifiable = True
model.settings.maxCo = 1
model.settings.maxDeltaT = 0.1

model.add_solver(thSolver)


#==============================================================================*
# Run

# Export to OpenFOAM format
ffn.allclean()
model.export_to_openfoam()

# Preprocessing
ffn.run_preprocessing(model=model)

model.plot_mesh(region=thMesh, show_edges=True)
model.plot_mesh(region=thMesh, show_edges=True, normal='y')
model.plot_mesh(region=thMesh, show_edges=True, normal='z')

model.plot_boundary(region=thMesh, boundaryName="pipeWall")
model.plot_boundary(region=thMesh, boundaryName="defaultFaces")

# Run the simulation
ffn.run(model=model)


#==============================================================================*
# Post-processing

model.plot_residuals(
    parameters=['p_rgh', 'h'],
    title="Thermal-hydraulics"
)
model.plot_mesh(
    region=thMesh,
    time=model.settings.endTime,
    fieldName='T',
    cmap="RdBu_r"
)
model.plot_mesh(
    region=thMesh,
    time=model.settings.endTime,
    fieldName='U',
    cmap="RdBu_r"
)
model.plot_mesh(
    region=thMesh,
    time=model.settings.endTime,
    fieldName='p_rgh',
    cmap="RdBu_r"
)


model.plot_animation(
    region=thMesh,
    fieldName="T",
    cmap="RdBu_r",
    unit='K',
    fps=10
)
model.plot_animation(
    region=thMesh,
    fieldName="U",
    cmap="RdBu_r",
    unit='m/s',
    fps=10
)

#==============================================================================*
