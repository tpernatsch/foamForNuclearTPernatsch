import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

nMesh.create_half_sphere(
    "center",
    radius=1,
    nCenter=4, nBorder=4,
    isAddAllBC=True
)
nMesh.create_hollow_half_sphere(
    "ringSphere1",
    innerRadius=1, outerRadius=2,
    nr=2, nt=4,
    isAddAllBC=True
)
nMesh.create_hollow_half_sphere(
    "ringSphere2",
    innerRadius=2, outerRadius=3,
    nr=1, nt=4,
    isAddAllBC=True
)
nMesh.create_cylinder_along_z(
    name="core",
    radius=1,
    lowZ=0, highZ=1,
    nx=4, ny=4, nz=4,
    isAddAllBC=True
)
nMesh.create_ring_along_z(
    name="ring",
    innerRadius=2,
    outerRadius=3,
    lowZ=0, highZ=1,
    nr=1, nt=4, nz=3,
    isAddAllBC=True
)

nMesh.add_merge_patch_pairs()
nMesh.merge_patches_with_name(name='top', includeFacename=['Top_'])
nMesh.merge_patches_with_name(name='wall', includeFacename=['Wall'])

solver = ffn.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
for boundaryName in ['top', 'wall']:
    model.plot_boundary(region=nMesh, boundaryName=boundaryName, show_edges=True)
