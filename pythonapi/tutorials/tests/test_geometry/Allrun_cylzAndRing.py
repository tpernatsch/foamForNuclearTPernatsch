import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

nMesh.create_cylinder_along_z(
    "cylinder",
    radius=1,
    lowZ=0, highZ=1,
    nx=4, ny=4, nz=2,
    isAddAllBC=True
)
nMesh.create_ring_along_z(
    "ring",
    innerRadius=1, outerRadius=2, lowZ=0, highZ=1,
    nr=3, nt=4, nz=2,
    isAddAllBC=True
)

nMesh.add_merge_patch_pairs()
nMesh.merge_patches_with_name(name='outerClad', includeFacename=['OuterWall'])
nMesh.merge_patches_with_name(name='topFuel', includeFacename=['cylinderTop_'])
nMesh.merge_patches_with_name(name='bottomFuel', includeFacename=['cylinderBottom_'])
nMesh.merge_patches_with_name(name='topClad', includeFacename=['ringTop_'])
nMesh.merge_patches_with_name(name='bottomClad', includeFacename=['ringBottom_'])

solver = ffn.solvers.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.case.Case()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
for boundaryName in ['topFuel', 'bottomFuel', 'outerClad', 'topClad', 'bottomClad']:
    model.plot_boundary(region=nMesh, boundaryName=boundaryName, show_edges=True)
