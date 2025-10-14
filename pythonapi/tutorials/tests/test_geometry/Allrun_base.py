import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

"""
nMesh.createCube(
    name='block',
    lowX=0, highX=1,
    lowY=0, highY=1,
    lowZ=0, highZ=1,
    nx=2, ny=3, nz=4
)

nMesh.createWedge(
    name='block',
    innerRadius=0,
    outerRadius=0.5,
    lowZ=0, highZ=1,
    wedgeAngle=15,
    nr=4,
    nz=5
)
"""

# nMesh.createCylinderAlongZ(
#     name='block',
#     radius=0.5,
#     lowZ=0, highZ=1,
#     nx=3, ny=3, nz=4,
#     x=1, y=-1
# )

# nMesh.createRingAlongZ(
#     name='block',
#     innerRadius=0.2, outerRadius=0.5,
#     lowZ=0, highZ=1,
#     nr=3, nt=10, nz=4,
#     x=1, y=-1
# )

# nMesh.createTriangularChannel(
#     name='block',
#     lowZ=0, highZ=1,
#     pitch=1,
#     radius=0.3,
#     nx=3, ny=3, nz=4
# )

# nMesh.createSphere(
#     name="sphere",
#     radius=1,
#     nCenter=10,
#     nBorder=2,
#     isAddBoundaryConditions=True
# )
# nMesh.mergePatchesWithName(name="wall", regex=".*Wall.*", patchType="wall")

# nMesh.createHollowHalfSphere(
#     name="sph1",
#     innerRadius=0.5, outerRadius=0.8,
#     nr=3,
#     nt=10,
#     isAddBoundaryConditions=True
# )

# nMesh.createHalfSphere(
#     name="halfSphere",
#     radius=1,
#     z=-1,
#     nCenter=10,
#     nBorder=10,
#     isAddBoundaryConditions=True
# )

nMesh.createSphere1D(
    name='core',
    innerRadius=0, outerRadius=1,
    opening=5,
    nr=10,
    isAddBoundaryConditions=True
)
# nMesh.createSphere1D(
#     name='reflector',
#     innerRadius=1, outerRadius=1.2,
#     opening=5,
#     nr=5,
#     isAddBoundaryConditions=True
# )
# nMesh.addMergePatchPairs(includeFacename=['Wall'], excludeFacename=['Wedge'])


# nMesh.createHalfSphere(
#     name="halfSphere",
#     radius=1,
#     z=-1,
#     nCenter=10,
#     nBorder=10,
#     isAddBoundaryConditions=True
# )
# nMesh.createRingAlongZ(
#     name="ring",
#     innerRadius=0.8,
#     outerRadius=1,
#     lowZ=-1, highZ=1,
#     nr=3, nt=10, nz=4,
#     isAddBoundaryConditions=True
# )
# nMesh.createCylinderAlongZ(
#     name="cylz",
#     radius=0.5, lowZ=-1, highZ=1,
#     nx=10, ny=10, nz=4,
#     isAddBoundaryConditions=True
# )
# nMesh.mergeBoundaryFaces(
#     includeFacename=["ringBottom", "cylzBottom"],
#     newBoundaryName="mergedRingsBottom",
#     newBoundaryType="wall"
# )
# nMesh.mergePatchPairs.append(("mergedRingsBottom", "halfSphereWallTop_0"))


solver = ffn.NeutronicsSolver(
    region=nMesh.region,
    mesh=nMesh,
    solver='diffusionNeutronics'
)


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
try:
    model.plot_boundary(region=nMesh, boundaryName="wall", show_edges=True)
except:
    pass
