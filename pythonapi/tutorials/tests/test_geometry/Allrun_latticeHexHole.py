import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

pitch = 0.01

# nMesh.createHexagonPrismWithHoleAlongZ(
#     "hexHole",
#     zmin=0, zmax=0.01,
#     pitch=pitch,
#     radius=0.003,
#     nr=4, nt=4, nz=2,
#     isAddBoundaryConditions=True
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

nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createHexagonPrismWithHoleAlongZ(
        "hexHole",
        zmin=0, zmax=0.01,
        pitch=pitch,
        radius=0.003,
        x=x, y=y,
        nr=4, nt=4, nz=2,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType='hexagon',
    pitch=pitch,
    nx=latticeNX,
    ny=latticeNY,
    elementsToPlace='F'
)
nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createHexagonPrismWithHoleAlongZ(
        "hexHole",
        zmin=0, zmax=0.01,
        pitch=pitch,
        radius=0.001,
        x=x, y=y,
        nr=7, nt=4, nz=2,
        isAddBoundaryConditions=True
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

nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createEdgeHexagonPrismAlongZ(
        "hexEdge",
        zmin=0, zmax=0.01,
        pitch=pitch,
        distanceCenterToEdge=0.002,
        edgeFaceOrientation=angle(x, y),
        x=x, y=y,
        nr=3, nt=4, nz=2,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType='hexagon',
    pitch=pitch,
    nx=latticeNX,
    ny=latticeNY,
    elementsToPlace='E'
)

nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createCornerHexagonPrismAlongZ(
        "hexCorner",
        zmin=0, zmax=0.01,
        pitch=pitch,
        distanceCenterToEdge=0.002,
        edgeFaceOrientation=np.arctan2(y, x)-np.pi/6,
        x=x, y=y,
        nr=3, nt=4, nz=2,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType='hexagon',
    pitch=pitch,
    nx=latticeNX,
    ny=latticeNY,
    elementsToPlace='C'
)

nMesh.addMergePatchPairs()
# nMesh.mergePatchesWithName(name='outerClad', includeFacename=['OuterWall'])
nMesh.mergePatchesWithName(name='top', includeFacename=['Top'])
nMesh.mergePatchesWithName(name='bottom', includeFacename=['Bottom'])
nMesh.mergePatchesWithName(name='wall', includeFacename=['Wall'])

solver = ffn.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.Model()
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
