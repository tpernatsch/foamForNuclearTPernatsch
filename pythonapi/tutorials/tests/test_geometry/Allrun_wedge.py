import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMeshWedge(appertureAngle=20, region='neutroMesh')

block1 = nMesh.createWedge(
    "block1",
    innerRadius=0,
    outerRadius=1,
    lowZ=0,
    highZ=1,
    nr=2,
    nz=3
)
block21 = nMesh.extrudeRight([block1], "block2", dr=0.5, nr=2)

block22 = nMesh.extrudeTop([block21], "block2", dz=1, nz=3)
block23 = nMesh.extrudeTop([block22], "block2", dz=1, nz=3)

block13 = nMesh.extrudeLeft([block23], "block13", dr=1, nr=3)

(block31, block32, block33) = nMesh.extrudeRight([block21, block22, block23], "block3", dr=1, nr=2)


mesh.BlockMesh.extrudeRight(nMesh, [block33], 'new', 1.3, 10)


solver = ffn.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
model.plot_mesh(region=nMesh, show_edges=True, normal='y')
model.plot_mesh(region=nMesh, show_edges=True, normal='z')
