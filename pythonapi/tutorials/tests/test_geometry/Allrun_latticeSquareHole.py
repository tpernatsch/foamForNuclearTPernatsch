import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

# nMesh.createCubeWithHoleAlongZ(
#     "negCylz",
#     lowX=0, highX=0.01,
#     lowY=0, highY=0.01,
#     lowZ=0, highZ=0.01,
#     radius=0.003,
#     nx=4, ny=4, nz=2, nt=4,
#     isAddBoundaryConditions=True
# )

pitch = 0.01
lattice = """F F
F F
"""

nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createCubeWithHoleAlongZ(
        "squareHole",
        lowX=x-pitch/2, highX=x+pitch/2,
        lowY=y-pitch/2, highY=y+pitch/2,
        lowZ=0, highZ=0.01,
        radius=0.003,
        nx=4, ny=4, nz=2, nt=4,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType='square',
    pitch=pitch,
    nx=2,
    ny=2
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
model.plot_boundary(region=nMesh, boundaryName='top', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='bottom', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='wall', show_edges=True)
