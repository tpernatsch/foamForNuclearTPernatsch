import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

pitch = 0.01

# nMesh.create_hexagon_prism_with_hole_along_z(
#     "hexHole",
#     zmin=0, zmax=0.01,
#     pitch=pitch,
#     radius=0.003,
#     nr=4, nt=4, nz=2,
#     isAddAllBC=True
# )

lattice = """0 0 0 C E E C
 0 0 E W F W E
  0 E F F F F E
   C W F F F W C
    E F F F F E 0
     E W F W E 0 0
      C E E C 0 0 0"""
latticeNX = len(lattice.split("\n"))
latticeNY = latticeNX

nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_hexagon_prism_with_hole_along_z(
        "hexHole",
        zmin=0, zmax=0.01,
        pitch=pitch,
        radius=0.003,
        x=x, y=y,
        nr=4, nt=4, nz=2,
        isAddAllBC=True
    ),
    lattice=lattice,
    latticeType='hexagon',
    pitch=pitch,
    nx=latticeNX,
    ny=latticeNY,
    elementsToPlace='F'
)
nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_hexagon_prism_with_hole_along_z(
        "hexHole",
        zmin=0, zmax=0.01,
        pitch=pitch,
        radius=0.001,
        x=x, y=y,
        nr=7, nt=4, nz=2,
        isAddAllBC=True
    ),
    lattice=lattice,
    latticeType='hexagon',
    pitch=pitch,
    nx=latticeNX,
    ny=latticeNY,
    elementsToPlace='W'
)

def angle(x, y):
    theta = np.arctan2(y, x)  # Angle in radians, range [-π, π]
    if theta < 0:
        theta += 2 * np.pi  # Normalize to [0, 2π)

    quantized = round(theta / (np.pi / 3)) * (np.pi / 3)
    return quantized

nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_edge_hexagon_prism_along_z(
        "hexEdge",
        zmin=0, zmax=0.01,
        pitch=pitch,
        distanceCenterToEdge=0.002,
        edgeFaceOrientation=angle(x, y),
        x=x, y=y,
        nr=3, nt=4, nz=2,
        isAddAllBC=True
    ),
    lattice=lattice,
    latticeType='hexagon',
    pitch=pitch,
    nx=latticeNX,
    ny=latticeNY,
    elementsToPlace='E'
)

nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_corner_hexagon_prism_along_z(
        "hexCorner",
        zmin=0, zmax=0.01,
        pitch=pitch,
        distanceCenterToEdge=0.002,
        edgeFaceOrientation=np.arctan2(y, x)-np.pi/6,
        x=x, y=y,
        nr=3, nt=4, nz=2,
        isAddAllBC=True
    ),
    lattice=lattice,
    latticeType='hexagon',
    pitch=pitch,
    nx=latticeNX,
    ny=latticeNY,
    elementsToPlace='C'
)

nMesh.add_merge_patch_pairs()
# nMesh.merge_patches_with_name(name='outerClad', includeFacename=['OuterWall'])
nMesh.merge_patches_with_name(name='top', includeFacename=['Top'])
nMesh.merge_patches_with_name(name='bottom', includeFacename=['Bottom'])
nMesh.merge_patches_with_name(name='wall', includeFacename=['Wall'])

solver = ffn.solvers.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.case.Case()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='top', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='bottom', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='wall', show_edges=True)
