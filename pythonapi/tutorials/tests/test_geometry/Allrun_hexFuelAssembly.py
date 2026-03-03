import numpy as np
import foamForNuclear as ffn
import foamForNuclear.mesh as mesh


nMesh = mesh.BlockMesh(region='neutroMesh')


lattice = """
0 0 0 C E E C
 0 0 E W F W E
  0 E F F F F E
   C W F F F W C
    E F F F F E 0
     E W F W E 0 0
      C E E C 0 0 0
"""
latticeNXY = len(lattice.strip().split("\n"))
# latticeNY = latticeNX

pitch = 0.01
zmin, zmax = 0, 0.01
nt = 4
nz = 2

funcPin = {
    'F': lambda name, x, y: nMesh.create_hexagon_prism_with_hole_along_z(
        name,
        zmin=zmin, zmax=zmax,
        pitch=pitch,
        radius=0.003,
        x=x, y=y,
        nr=3, nt=nt, nz=nz,
        isAddAllBC=True
    ),
    'W': lambda name, x, y: nMesh.create_hexagon_prism_with_hole_along_z(
        name,
        zmin=zmin, zmax=zmax,
        pitch=pitch,
        radius=0.001,
        x=x, y=y,
        nr=5, nt=nt, nz=nz,
        isAddAllBC=True
    )
}

assemblyPitch = 0.025*2


nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.hexagonal_lattice_assembly(
        pitch, lattice, latticeNXY, zmin, zmax, assemblyPitch,
        nrEdge=2, ntEdge=nt, nzEdge=nz,
        xAssembly=x, yAssembly=y,
        funcElementGeneratorPin=funcPin
    ),
    lattice="""0 0 0
    0 F F
    0 F 0""",
    nx=3,
    ny=3,
    pitch=assemblyPitch,
    latticeType='hexagon',
    flatToFlatDirection='x',
    elementsToPlace='F'
)


nMesh.add_merge_patch_pairs()
# nMesh.merge_patches_with_name(name='outerClad', includeFacename=['OuterWall'])
nMesh.merge_patches_with_name(name='top', includeFacename=['Top'])
nMesh.merge_patches_with_name(name='bottom', includeFacename=['Bottom'])
nMesh.merge_patches_with_name(name='wall', includeFacename=['Wall'])

solver = ffn.solvers.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.case.Case()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=False)
# model.plot_boundary(region=nMesh, boundaryName='defaultFaces', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='top', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='bottom', show_edges=True)
# model.plot_boundary(region=nMesh, boundaryName='wall', show_edges=True)
