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

lattice = """0 0 0 0 0 0 0 R R R R R R R R
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
              R R R R R R R R 0 0 0 0 0 0 0"""

latticeNX = len([line for line in lattice.split("\n") if line != ""])
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

nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createHexagonPrismAlongZ(
        name="innerCore",
        zmin=-coreHeight/2,
        zmax=coreHeight/2,
        pitch=flatToFlat,
        x=x, y=y,
        nr=nr,
        nt=nt,
        nz=coreNodes,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType="hexagon",
    nx=latticeNY, ny=latticeNY, pitch=pitch,
    elementsToPlace=['I']
)
nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createHexagonPrismAlongZ(
        name="outerCore",
        zmin=-coreHeight/2,
        zmax=coreHeight/2,
        pitch=flatToFlat,
        x=x, y=y,
        nr=nr,
        nt=nt,
        nz=coreNodes,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType="hexagon",
    nx=latticeNY, ny=latticeNY, pitch=pitch,
    elementsToPlace=['O']
)
nMesh.latticePlacement(
    funcElementGenerator=lambda name, x, y: nMesh.createHexagonPrismAlongZ(
        name="upperGasPlenum",
        zmin=coreHeight/2,
        zmax=coreHeight/2+upperFGPHeight,
        pitch=flatToFlat,
        x=x, y=y,
        nr=nr,
        nt=nt,
        nz=upperFGPNodes,
        isAddBoundaryConditions=True
    ),
    lattice=lattice,
    latticeType="hexagon",
    nx=latticeNY, ny=latticeNY, pitch=pitch,
    elementsToPlace=['I', 'O']
)

nMesh.addMergePatchPairs(
    includeFacename=['Top', 'Bottom'],
    excludeFacename=['Wall']
)
if (isFluid):
    bafflesFaces = nMesh.addBaffles(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom']
    )
else:
    nMesh.addMergePatchPairs(
        includeFacename=['Wall'],
        excludeFacename=['Top', 'Bottom'],
    )

externalWalls = nMesh.getStandaloneFaces(
    includeFacename=['Wall'],
    excludeFacename=['Top', 'Bottom'],
)

nMesh.mergePatchesWithName(
    name="wall",
    includeFacename=[face.name for face in externalWalls],
    patchType="wall"
)
nMesh.mergePatchesWithName(name="bottom", includeFacename=["Bottom_"])
nMesh.mergePatchesWithName(name="top", includeFacename=["Top_"])


solver = ffn.NeutronicsSolver(region=nMesh.region, mesh=nMesh, solver='diffusionNeutronics')


model = ffn.Model()
model.settings.application = 'dummy'

model.solvers.append(solver)

ffn.allclean()
model.export_to_openfoam()

ffn.run_preprocessing(model=model)

model.plot_mesh(region=nMesh, show_edges=True)

model.plot_boundary(region=nMesh, boundaryName="wall", show_edges=True)
