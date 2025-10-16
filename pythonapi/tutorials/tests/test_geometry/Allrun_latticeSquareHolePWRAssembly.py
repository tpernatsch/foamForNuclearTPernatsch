import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')


pitch = 0.0126

lattice = """F F F F F F F F F F F F F F F F F
F F F F F F F F F F F F F F F F F
F F G F G C F G C G F C G F G F F
F F F C F F F F F F F F F C F F F
F F G F F F G F F F G F F F G F F
F F C F F C F F C F F C F F C F F
F F F F G F F F G F F F G F F F F
F F G F F F F F F F F F F F G F F
F F C F F C G F C F G C F F C F F
F F G F F F F F F F F F F F G F F
F F F F G F F F G F F F G F F F F
F F C F F C F F C F F C F F C F F
F F G F F F G F F F G F F F G F F
F F F C F F F F F F F F F C F F F
F F G F G C F G C G F C G F G F F
F F F F F F F F F F F F F F F F F
F F F F F F F F F F F F F F F F F"""

# lattice = """F F F
# F F F
# F F C"""

nxLat = len(lattice.split('\n'))
nyLat = nxLat

lowZ = 0
highZ = 0.1

nx = 4
ny = nx
nz = 3
nt = 4

claddingOR  =  9.500e-3 / 2
guideTudeOR = 12.243e-3 / 2

nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createCubeWithHoleAlongZ(
        "squareHole",
        lowX=x-pitch/2, highX=x+pitch/2,
        lowY=y-pitch/2, highY=y+pitch/2,
        lowZ=lowZ, highZ=highZ,
        radius=claddingOR,
        nx=nx, ny=ny, nz=nz, nt=nt,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType='square',
    pitch=pitch,
    nx=nxLat,
    ny=nyLat,
    elementsToPlace=['F', 'G']
)
nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createCubeWithHoleAlongZ(
        "guideTube",
        lowX=x-pitch/2, highX=x+pitch/2,
        lowY=y-pitch/2, highY=y+pitch/2,
        lowZ=lowZ, highZ=highZ,
        radius=guideTudeOR,
        nx=nx, ny=ny, nz=nz, nt=nt,
        squareEdgeToHoleCenter=claddingOR / np.sqrt(2),
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType='square',
    pitch=pitch,
    nx=nxLat,
    ny=nyLat,
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

model.plot_mesh(region=nMesh, show_edges=False)
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='top', show_edges=True, normal='z')
model.plot_boundary(region=nMesh, boundaryName='bottom', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='wall', show_edges=True)
