import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')

nMesh.createCylinderAlongZ(
    "cylinder",
    radius=1,
    lowZ=0, highZ=1,
    nx=4, ny=4, nz=2,
    isAddBoundaryConditions=True
)
nMesh.createRingAlongZ(
    "ring",
    innerRadius=1, outerRadius=2, lowZ=0, highZ=1,
    nr=3, nt=4, nz=2,
    isAddBoundaryConditions=True
)

nMesh.addMergePatchPairs()
nMesh.mergePatchesWithName(name='outerClad', includeFacename=['OuterWall'])
nMesh.mergePatchesWithName(name='topFuel', includeFacename=['cylinderTop_'])
nMesh.mergePatchesWithName(name='bottomFuel', includeFacename=['cylinderBottom_'])
nMesh.mergePatchesWithName(name='topClad', includeFacename=['ringTop_'])
nMesh.mergePatchesWithName(name='bottomClad', includeFacename=['ringBottom_'])

solver = ffn.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
for boundaryName in ['topFuel', 'bottomFuel', 'outerClad', 'topClad', 'bottomClad']:
    model.plot_boundary(region=nMesh, boundaryName=boundaryName, show_edges=True)
