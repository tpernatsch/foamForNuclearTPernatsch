import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMeshWedge(appertureAngle=20, region='neutroMesh')

block1 = nMesh.create_wedge(
    "block1",
    innerRadius=0,
    outerRadius=1,
    lowZ=0,
    highZ=1,
    nr=2,
    nz=3
)
block21 = nMesh.extrude_right([block1], "block2", dr=0.5, nr=2)

block22 = nMesh.extrude_top([block21], "block2", dz=1, nz=3)
block23 = nMesh.extrude_top([block22], "block2", dz=1, nz=3)

block13 = nMesh.extrude_left([block23], "block13", dr=1, nr=3)

(block31, block32, block33) = nMesh.extrude_right([block21, block22, block23], "block3", dr=1, nr=2)


mesh.BlockMesh.extrude_right(nMesh, [block33], 'new', 1.3, 10)


solver = ffn.solvers.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.case.Case()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)
model.plot_mesh(region=nMesh, show_edges=True, normal='y')
model.plot_mesh(region=nMesh, show_edges=True, normal='z')
