"""
mesh.py — Computational mesh for the 3D Small ESFR core (fluid region)
=======================================================================

The thermal-hydraulics solver works on the coolant channel network inside
the reactor core.  The mesh represents the sodium flow path from the inlet
at the bottom (diagrid plenum) to the outlet at the top.

Hexagonal assembly positions (same lattice as the neutronics mesh):
    I = inner core fuel assembly   O = outer core fuel assembly
    C = control rod assembly       R = radial reflector assembly

BAFFLES: The fluid mesh includes internal "baffles" — thin walls between
adjacent assemblies that force the sodium to flow upward through each
assembly independently (no cross-flow).  These are added automatically
when generating the mesh with isFluidMesh=True.

TWO MESH OPTIONS:
  Option A — load pre-generated polyMesh files (fast, recommended)
  Option B — generate mesh with blockMesh (~5-10 min)

Set USE_PREMADE_MESH below.
"""

import numpy as np
import foamForNuclear.mesh as mesh

# ---------------------------------------------------------------------------
USE_PREMADE_MESH = True
# ---------------------------------------------------------------------------

LATTICE = """0 0 0 0 0 0 0 R R R R R R R R
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

if USE_PREMADE_MESH:
    thMesh = mesh.PolyMesh(region="fluidRegion", srcpath="../meshes/polyMeshFluid")

else:
    EDGE             = 0.12171
    GAP              = 0.0
    CORE_HEIGHT      = 1.00
    LOW_REFL_HEIGHT  = 0.30
    LOW_FGP_HEIGHT   = 0.36
    DIAGRID_HEIGHT   = 0.30
    UPPER_FGP_HEIGHT = 0.11
    UPPER_REFL_HEIGHT= 0.38
    CORE_NODES       = 10
    LOW_REFL_NODES   = 3
    LOW_FGP_NODES    = 3
    DIAGRID_NODES    = 3
    UPPER_FGP_NODES  = 2
    UPPER_REFL_NODES = 3

    flat_to_flat = EDGE * np.sqrt(3)
    pitch        = flat_to_flat + GAP
    lattice_n    = len(LATTICE.split("\n"))
    nr, nt       = 1, 1

    thMesh = mesh.BlockMesh(region="fluidRegion")

    def _add_zone(name, zmin, zmax, nz, elements):
        thMesh.lattice_placement(
            funcElementGenerator=lambda n, x, y: thMesh.create_hexagon_prism_along_z(
                name=name, zmin=zmin, zmax=zmax,
                pitch=flat_to_flat, x=x, y=y,
                nr=nr, nt=nt, nz=nz, isAddAllBC=True
            ),
            lattice=LATTICE, latticeType="hexagon",
            nx=lattice_n, ny=lattice_n, pitch=pitch,
            elementsToPlace=elements,
        )

    H  = CORE_HEIGHT; UF = UPPER_FGP_HEIGHT; UR = UPPER_REFL_HEIGHT
    LF = LOW_FGP_HEIGHT; LR = LOW_REFL_HEIGHT; DG = DIAGRID_HEIGHT

    _add_zone("innerCore",      -H/2, H/2, CORE_NODES, ['I'])
    _add_zone("outerCore",      -H/2, H/2, CORE_NODES, ['O'])
    _add_zone("upperGasPlenum",  H/2,  H/2+UF,     UPPER_FGP_NODES,  ['I','O'])
    _add_zone("upperReflector",  H/2+UF, H/2+UF+UR, UPPER_REFL_NODES, ['I','O'])
    _add_zone("lowerGasPlenum", -H/2-LF, -H/2,      LOW_FGP_NODES,    ['I','O'])
    _add_zone("lowerReflector", -H/2-LF-LR, -H/2-LF, LOW_REFL_NODES,  ['I','O'])
    for zmin, zmax, nz in [
        (-H/2, H/2, CORE_NODES), (H/2, H/2+UF, UPPER_FGP_NODES),
        (H/2+UF, H/2+UF+UR, UPPER_REFL_NODES), (-H/2-LF, -H/2, LOW_FGP_NODES),
        (-H/2-LF-LR, -H/2-LF, LOW_REFL_NODES)
    ]:
        _add_zone("radialReflector", zmin, zmax, nz, ['R'])
    _add_zone("diagrid", -H/2-LF-LR-DG, -H/2-LF-LR, DIAGRID_NODES, ['I','O','R','C'])
    for zmin, zmax, nz in [
        (-H/2-LF-LR, -H/2-LF, LOW_REFL_NODES),
        (-H/2-LF,    -H/2,     LOW_FGP_NODES),
        (-H/2,        H/2,     CORE_NODES),
    ]:
        _add_zone("follower", zmin, zmax, nz, ['C'])
    _add_zone("controlRod", H/2,    H/2+UF,    UPPER_FGP_NODES,  ['C'])
    _add_zone("controlRod", H/2+UF, H/2+UF+UR, UPPER_REFL_NODES, ['C'])

    # Connect top/bottom faces of adjacent axial zones
    thMesh.add_merge_patch_pairs(includeFacename=['Top','Bottom'], excludeFacename=['Wall'])
    # Fluid mesh: add baffles between assemblies (prevent cross-flow)
    thMesh.add_baffles(includeFacename=['Wall'], excludeFacename=['Top','Bottom'])

    ext_walls = thMesh.get_standalone_faces(
        includeFacename=['Wall'], excludeFacename=['Top','Bottom']
    )
    thMesh.merge_patches_with_name(
        name="wall", includeFacename=[f.name for f in ext_walls], patchType="wall"
    )
    thMesh.merge_patches_with_name(name="bottom", includeFacename=["diagridBottom_"])
    thMesh.merge_patches_with_name(
        name="top",
        includeFacename=["upperReflectorTop_", "controlRodTop_", "radialReflectorTop_"]
    )
