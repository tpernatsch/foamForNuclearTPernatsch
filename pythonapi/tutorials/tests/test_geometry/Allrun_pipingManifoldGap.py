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

core = thMesh.create_cylinder_along_z(
    name='core',
    radius=coreRadius,
    lowZ=0, highZ=coreLength,
    nx=3, ny=3, nz=10,
    isAddAllBC=False
)
bottomPlenum = thMesh.create_half_sphere(
    name='bottomPlenum',
    radius=vesselRadius,
    nCenter=7, nBorder=3,
    isAddAllBC=False
)
ring = thMesh.create_ring_along_z(
    name='ring',
    innerRadius=coreRadius+gap,
    outerRadius=vesselRadius,
    lowZ=0, highZ=coreLength,
    nr=2, nt=7, nz=4,
    isAddAllBC=False
)

inletManifold = thMesh.create_pipe_cylindrical_manifold_along_z(
    'inletManifold',
    nEntries=nLoops,
    innerRadius=coreRadius+gap,
    outerRadius=vesselRadius,
    equivalentHydraulicDiameter=equivalentHydraulicDiameter,
    lowZ=coreLength,
    nr=2, nt=3,
    isAddGap=True
)
outletManifold = thMesh.create_pipe_cylindrical_manifold_along_z(
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
    thMesh.extrude_normal(inletPipe, 'front', f'inletPipe{i}', length=legLength, n=3)

for i, outletPipe in enumerate(outletManifold[::2]):
    thMesh.extrude_normal(outletPipe, 'front', f'outletPipe{i}', length=legLength + vesselRadius-coreRadius, n=3)


# Add inlet and outlet of manifold
inletManifoldOutlet = ffn.Face("inletManifoldOutlet")
for block in inletManifold:
    inletManifoldOutlet.add_sub_face(block.bottomFace())

outletManifoldInlet = ffn.Face("outletManifoldInlet")
for block in outletManifold:
    outletManifoldInlet.add_sub_face(block.bottomFace())

ringTop = ffn.Face("ringTop")
for block in ring:
    ringTop.add_sub_face(block.topFace())

coreTop = ffn.Face("coreTop")
for block in core:
    coreTop.add_sub_face(block.topFace())

coreBottom = ffn.Face("coreBottom")
for block in core:
    coreBottom.add_sub_face(block.bottomFace())
for block in ring:
    coreBottom.add_sub_face(block.bottomFace())

bottomPlenumTop = ffn.Face("bottomPlenumTop")
for block in bottomPlenum[:-1]:
    bottomPlenumTop.add_sub_face(block.topFace())

thMesh.add_boundary(inletManifoldOutlet)
thMesh.add_boundary(outletManifoldInlet)
thMesh.add_boundary(ringTop)
thMesh.add_boundary(coreTop)
thMesh.add_boundary(coreBottom)
thMesh.add_boundary(bottomPlenumTop)

# Stich blocks
thMesh.merge_patch_pairs_by_name(coreTop.name, outletManifoldInlet.name)
thMesh.merge_patch_pairs_by_name(ringTop.name, inletManifoldOutlet.name)
thMesh.merge_patch_pairs_by_name(coreBottom.name, bottomPlenumTop.name)


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
