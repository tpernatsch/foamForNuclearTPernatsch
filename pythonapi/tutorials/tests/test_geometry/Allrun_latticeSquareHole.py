import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

# nMesh.create_cube_with_hole_along_z(
#     "negCylz",
#     lowX=0, highX=0.01,
#     lowY=0, highY=0.01,
#     lowZ=0, highZ=0.01,
#     radius=0.003,
#     nx=4, ny=4, nz=2, nt=4,
#     isAddAllBC=True
# )

pitch = 0.01
lattice = """
F F
F F
"""

nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_cube_with_hole_along_z(
        "squareHole",
        lowX=x-pitch/2, highX=x+pitch/2,
        lowY=y-pitch/2, highY=y+pitch/2,
        lowZ=0, highZ=0.01,
        radius=0.003,
        nx=4, ny=4, nz=2, nt=4,
        isAddAllBC=True
    ),
    lattice=lattice.strip(),
    latticeType='square',
    pitch=pitch,
    nx=2,
    ny=2
)

nMesh.merge_patches_with_name(name='outerClad', includeFacename=['WallHole'])
nMesh.merge_patches_with_name(name='top', includeFacename=['Top'])
nMesh.merge_patches_with_name(name='bottom', includeFacename=['Bottom'])
externalWalls = nMesh.get_standalone_faces(
    includeFacename=['Wall'],
    excludeFacename=['Top', 'Bottom'],
)
nMesh.merge_patches_with_name(
    name="wall",
    includeFacename=[face.name for face in externalWalls],
    patchType="wall",
    isStrict=True
)

nMesh.isMergeCoincidentPoints = True

solver = ffn.solvers.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.case.Case()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='top', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='bottom', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='wall', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='outerClad', show_edges=True)
