"""

"""
#==============================================================================*
# Imports

import numpy as np
import foamForNuclear as ffn
import foamForNuclear.boundaryConditions as bc
import foamForNuclear.mesh as mesh


#==============================================================================*
# meshHex

lattice = """
0 0 0 0 0 0 0 R R R R R R R R
 0 0 0 0 0 0 R O C O O O C O R
  0 0 0 0 0 R C O O O O O O C R
   0 0 0 0 R O O I I I I I O O R
    0 0 0 R O O I C I I C I O O R
     0 0 R O O I I I I I I I O O R
      0 R C O I I I I I I I I O C R
       R O O I C I I I I I C I O O R
        R C O I I I I I I I I O C R 0
         R O O I I I I I I I O O R 0 0
          R O O I C I I C I O O R 0 0 0
           R O O I I I I I O O R 0 0 0 0
            R C O O O O O O C R 0 0 0 0 0
             R O C O O O C O R 0 0 0 0 0 0
              R R R R R R R R 0 0 0 0 0 0 0
              """

latticeNX = len([line for line in lattice.strip().split("\n") if line != ""])
latticeNY = latticeNX
nr = 1
nt = 1
edge = 0.12171 # 0.11391
gap = 0
flatToFlat = edge * np.sqrt(3)
pitch = flatToFlat + gap
coreHeight = 1
coreNodes = 10
upperFGPHeight = 0.11
upperFGPNodes = 2

isFluid = False

nMesh = mesh.BlockMesh(region='neutroRegion')

nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_hexagon_prism_along_z(
        name="innerCore",
        zmin=-coreHeight/2,
        zmax=coreHeight/2,
        pitch=flatToFlat,
        x=x, y=y,
        nr=nr,
        nt=nt,
        nz=coreNodes,
        isAddLateralBC=True,
        isAddBottomBC=True,
        isAddAllBC=isFluid
    ),
    lattice=lattice,
    latticeType="hexagon",
    nx=latticeNY, ny=latticeNY, pitch=pitch,
    elementsToPlace=['I']
)
nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_hexagon_prism_along_z(
        name="outerCore",
        zmin=-coreHeight/2,
        zmax=coreHeight/2,
        pitch=flatToFlat,
        x=x, y=y,
        nr=nr,
        nt=nt,
        nz=coreNodes,
        isAddLateralBC=True,
        isAddBottomBC=True,
        isAddAllBC=isFluid
    ),
    lattice=lattice,
    latticeType="hexagon",
    nx=latticeNY, ny=latticeNY, pitch=pitch,
    elementsToPlace=['O']
)
nMesh.lattice_placement(
    funcElementGenerator=lambda name, x, y: nMesh.create_hexagon_prism_along_z(
        name="upperGasPlenum",
        zmin=coreHeight/2,
        zmax=coreHeight/2+upperFGPHeight,
        pitch=flatToFlat,
        x=x, y=y,
        nr=nr,
        nt=nt,
        nz=upperFGPNodes,
        isAddLateralBC=True,
        isAddTopBC=True,
        isAddAllBC=isFluid
    ),
    lattice=lattice,
    latticeType="hexagon",
    nx=latticeNY, ny=latticeNY, pitch=pitch,
    elementsToPlace=['I', 'O']
)

if (isFluid):
    # Use mergePatchPairs for fluid regions with baffles
    nMesh.add_merge_patch_pairs(
        includeFacename=['Top', 'Bottom'],
        excludeFacename=['Wall']
    )
    bafflesFaces = nMesh.add_baffles(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom']
    )
else:
    nMesh.isMergeCoincidentPoints = True

# Create BCs
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

nMesh.merge_patches_with_name(name="bottom", includeFacename=["Bottom_"])
nMesh.merge_patches_with_name(name="top", includeFacename=["Top_"])



solver = ffn.NeutronicsSolver(
    region=nMesh.region,
    mesh=nMesh,
    solver='diffusionNeutronics'
)


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)

# model.plot_boundary(region=nMesh, boundaryName="wall", show_edges=True)
