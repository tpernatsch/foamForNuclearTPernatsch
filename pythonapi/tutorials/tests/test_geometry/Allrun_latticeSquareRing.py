import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

pitch = 0.01

# lattice = """
# 0 F F F F 0
# F F F F F F
# F F F F F F
# F F F F F F
# F F F F F F
# 0 F F F F 0
# """
lattice = """
0 0 0 0 K K K K K K 0 0 0 0
0 0 K K K F F F F K K K 0 0
0 K K F F F F F F F F K K 0
0 K F F F F F F F F F F K 0
K K F F F F F F F F F F K K
K F F F F F F F F F F F F K
K F F F F F F F F F F F F K
K F F F F F F F F F F F F K
K F F F F F F F F F F F F K
K K F F F F F F F F F F K K
0 K F F F F F F F F F F K 0
0 K K F F F F F F F F K K 0
0 0 K K K F F F F K K K 0 0
0 0 0 0 K K K K K K 0 0 0 0
"""

nXY = len(lattice.strip('\n').split('\n'))

nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_cube_with_hole_along_z(
        "squareHole",
        lowX=x-pitch/2, highX=x+pitch/2,
        lowY=y-pitch/2, highY=y+pitch/2,
        lowZ=0, highZ=0.01,
        radius=0.003,
        nx=4, ny=4, nz=2, nt=4,
        isAddAllBC=True,
    ),
    lattice=lattice,
    latticeType='square',
    pitch=pitch,
    nx=nXY,
    ny=nXY,
    elementsToPlace='F'
)
nMesh.fill_lattice_ring_gap(
    name='gap',
    # ringRadius=0.039,
    ringRadius=0.065,
    nxLat=nXY,
    nyLat=nXY,
    pitch=pitch,
    latticeType='square',
    nxBlock=4,
    nyBlock=4,
    nzBlock=2,
    zmin=0,
    zmax=0.01,
    isAddAllBC=True
)

# nMesh.add_merge_patch_pairs()
nMesh.merge_patches_with_name(name='top', includeFacename=['Top'])
nMesh.merge_patches_with_name(name='bottom', includeFacename=['Bottom'])
nMesh.merge_patches_with_name(name='channelWall', includeFacename=['squareHoleWallHole'])
nMesh.merge_patches_with_name(name='ringWall', includeFacename=['gapWallOuter'])

nMesh.isMergeCoincidentPoints = True

solver = ffn.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
model.plot_mesh(region=nMesh, show_edges=True, normal='z')
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='top', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='bottom', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='wall', show_edges=True)
