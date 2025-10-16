#==============================================================================*
# Imports

import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


#==============================================================================*
# Mesh

thMesh = mesh.BlockMesh(region='fluidRegion')

nLoops = 3

vesselRadius = 2.75
coreRadius = 2.2
gap = 0.1

coreLength = 4.6
legLength = 3

elbowRadius = 0.5
equivalentHydraulicDiameter = 0.6

core = thMesh.createCylinderAlongZ(
    name='core',
    radius=coreRadius,
    lowZ=0, highZ=coreLength,
    nx=3, ny=3, nz=10,
    isAddBoundaryConditions=False
)
bottomPlenum = thMesh.createHalfSphere(
    name='bottomPlenum',
    radius=vesselRadius,
    nCenter=7, nBorder=3,
    isAddBoundaryConditions=False
)
ring = thMesh.createRingAlongZ(
    name='ring',
    innerRadius=coreRadius+gap,
    outerRadius=vesselRadius,
    lowZ=0, highZ=coreLength,
    nr=2, nt=7, nz=4,
    isAddBoundaryConditions=False
)

inletManifold = thMesh.createPipeCylindricalManifoldAlongZ(
    'inletManifold',
    nEntries=nLoops,
    innerRadius=coreRadius+gap,
    outerRadius=vesselRadius,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    lowZ=coreLength,
    nr=2, nt=3,
    isAddGap=True
)
outletManifold = thMesh.createPipeCylindricalManifoldAlongZ(
    'outletManifold',
    nEntries=nLoops,
    innerRadius=0.001,
    outerRadius=coreRadius,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    lowZ=coreLength,
    angleStart=360/(2*nLoops),
    nr=3, nt=5,
    isAddGap=False
)

# Add pipes
for i, inletPipe in enumerate(inletManifold[:nLoops]):
    thMesh.extrudeNormal(inletPipe, 'front', f'inletPipe{i}', length=legLength, n=3)

for i, outletPipe in enumerate(outletManifold[::2]):
    thMesh.extrudeNormal(outletPipe, 'front', f'outletPipe{i}', length=legLength + vesselRadius-coreRadius, n=3)


# Add inlet and outlet of manifold
inletManifoldOutlet = ffn.Face("inletManifoldOutlet")
for block in inletManifold:
    inletManifoldOutlet.addSubFace(block.bottomFace())

outletManifoldInlet = ffn.Face("outletManifoldInlet")
for block in outletManifold:
    outletManifoldInlet.addSubFace(block.bottomFace())

ringTop = ffn.Face("ringTop")
for block in ring:
    ringTop.addSubFace(block.topFace())

coreTop = ffn.Face("coreTop")
for block in core:
    coreTop.addSubFace(block.topFace())

coreBottom = ffn.Face("coreBottom")
for block in core:
    coreBottom.addSubFace(block.bottomFace())
for block in ring:
    coreBottom.addSubFace(block.bottomFace())

bottomPlenumTop = ffn.Face("bottomPlenumTop")
for block in bottomPlenum[:-1]:
    bottomPlenumTop.addSubFace(block.topFace())

thMesh.addBoundary(inletManifoldOutlet)
thMesh.addBoundary(outletManifoldInlet)
thMesh.addBoundary(ringTop)
thMesh.addBoundary(coreTop)
thMesh.addBoundary(coreBottom)
thMesh.addBoundary(bottomPlenumTop)

# Stich blocks
thMesh.mergePatchPairsByName(coreTop.name, outletManifoldInlet.name)
thMesh.mergePatchPairsByName(ringTop.name, inletManifoldOutlet.name)
thMesh.mergePatchPairsByName(coreBottom.name, bottomPlenumTop.name)


#==============================================================================*
# Solver

thSolver = ffn.ThermalHydraulicsSolver(
    region=thMesh.region,
    mesh=thMesh,
    solver='onePhase',
    isSetFvSchemesToDefault=False
)


#==============================================================================*
# Model

model = ffn.Model()
model.settings.application = "dummy"

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


#==============================================================================*
