import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

radius = 1
pitch = 2.5
lowZ = 0
highZ = 1
nx = 3
ny = 4
nz = 2

# nMesh.create_triangular_channel(
#     name='fluid',
#     lowZ=lowZ, highZ=highZ,
#     pitch=pitch, radius=radius,
#     thetaZ=60*np.pi/180,
#     x=0, y=0,
#     nx=nx, ny=ny, nz=nz,
#     isAddAllBC=True
# )


lattice = """0 0 F F F
 0 F F F F
  F F F F F
   F F F F 0
    F F F 0 0"""
latticeNXY = len(lattice.split("\n"))
# latticeNY = latticeNX

nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_triangular_channel(
        name='fluid',
        lowZ=lowZ, highZ=highZ,
        pitch=pitch, radius=radius,
        thetaZ=-30*np.pi/180,
        x=x, y=y,
        nx=nx, ny=ny, nz=nz,
        isAddAllBC=True
    ),
    lattice=lattice,
    nx=latticeNXY,
    ny=latticeNXY,
    pitch=pitch,
    latticeType='hexagon',
    elementsToPlace='F'
)
nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_triangular_channel(
        name='fluid',
        lowZ=lowZ, highZ=highZ,
        pitch=pitch, radius=radius,
        thetaZ=30*np.pi/180,
        x=x, y=y,
        nx=nx, ny=ny, nz=nz,
        isAddAllBC=True
    ),
    lattice=lattice,
    nx=latticeNXY,
    ny=latticeNXY,
    pitch=pitch,
    x=pitch/np.sqrt(3),
    latticeType='hexagon',
    elementsToPlace='F'
)



nMesh.add_merge_patch_pairs()
nMesh.merge_patches_with_name(name='outerClad', includeFacename=['WallHole'])
nMesh.merge_patches_with_name(name='fluidPatch', includeFacename=['WallFluid'])
nMesh.merge_patches_with_name(name='outlet', includeFacename=['Top_'])
nMesh.merge_patches_with_name(name='inlet', includeFacename=['Bottom_'])

solver = ffn.solvers.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.case.Case()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
model.plot_mesh(region=nMesh, show_edges=True, normal='z')
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='inlet', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='outlet', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='outerClad', show_edges=True)
model.plot_boundary(region=nMesh, boundaryName='fluidPatch', show_edges=True)
