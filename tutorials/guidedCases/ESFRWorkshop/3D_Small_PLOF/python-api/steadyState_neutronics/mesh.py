"""
mesh.py — Computational mesh for the 3D Small ESFR core (neutronics region)
============================================================================

The reactor core is a hexagonal assembly lattice. Each position in the lattice
is one assembly type:

    I = inner core fuel assembly
    O = outer core fuel assembly
    C = control rod assembly (absorber + follower)
    R = radial reflector assembly
    0 = empty (outside core boundary)

The mesh is 3D: axially it spans the active core, lower/upper gas plenums,
axial reflectors, and the diagrid support structure below.

TWO MESH OPTIONS are available:
  Option A — Pre-generated mesh (recommended for the workshop):
             Loads existing mesh files from ../meshes/. Fast, no extra tools needed.
  Option B — Generate mesh from scratch:
             Runs the blockMesh algorithm to build the hexagonal prism mesh.
             Takes ~5-10 minutes. Useful if you want to change the geometry.

Set USE_PREMADE_MESH below to choose.
"""

import numpy as np
import foamForNuclear.mesh as mesh

# ---------------------------------------------------------------------------
USE_PREMADE_MESH = True   # True = Option A (fast), False = Option B (generate)
# ---------------------------------------------------------------------------

# Reactor lattice map (top view of the core)
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
    # -----------------------------------------------------------------------
    # Option A: load pre-generated polyMesh files
    # These were generated once with Option B and stored so you don't have to
    # wait each time you run the script.
    # -----------------------------------------------------------------------
    nMesh = mesh.PolyMesh(region="neutroRegion", srcpath="../meshes/polyMeshNeutro")

else:
    # -----------------------------------------------------------------------
    # Option B: generate the hexagonal prism mesh with blockMesh
    # Core geometry (all dimensions in metres)
    # -----------------------------------------------------------------------
    EDGE             = 0.12171  # assembly hex half-edge length [m]
    GAP              = 0.0      # inter-assembly gap [m]

    # Axial region heights [m]
    CORE_HEIGHT      = 1.00
    LOW_REFL_HEIGHT  = 0.30
    LOW_FGP_HEIGHT   = 0.36     # lower fission-gas plenum
    DIAGRID_HEIGHT   = 0.30
    UPPER_FGP_HEIGHT = 0.11
    UPPER_REFL_HEIGHT= 0.38

    # Number of axial mesh layers per region
    CORE_NODES       = 10
    LOW_REFL_NODES   = 3
    LOW_FGP_NODES    = 3
    DIAGRID_NODES    = 3
    UPPER_FGP_NODES  = 2
    UPPER_REFL_NODES = 3

    flat_to_flat = EDGE * np.sqrt(3)
    pitch        = flat_to_flat + GAP
    lattice_n    = len(LATTICE.split("\n"))
    nr, nt       = 1, 1   # no radial/tangential refinement

    nMesh = mesh.BlockMesh(region="neutroRegion")

    def _add_zone(name, zmin, zmax, nz, elements):
        """Place hexagonal prism blocks for a given zone type and axial range."""
        nMesh.lattice_placement(
            funcElementGenerator=lambda n, x, y: nMesh.create_hexagon_prism_along_z(
                name=name, zmin=zmin, zmax=zmax,
                pitch=flat_to_flat, x=x, y=y,
                nr=nr, nt=nt, nz=nz, isAddAllBC=True
            ),
            lattice=LATTICE, latticeType="hexagon",
            nx=lattice_n, ny=lattice_n, pitch=pitch,
            elementsToPlace=elements,
        )

    H  = CORE_HEIGHT
    UF = UPPER_FGP_HEIGHT
    UR = UPPER_REFL_HEIGHT
    LF = LOW_FGP_HEIGHT
    LR = LOW_REFL_HEIGHT
    DG = DIAGRID_HEIGHT

    # Active core zones
    _add_zone("innerCore",      -H/2,          H/2,           CORE_NODES,       ['I'])
    _add_zone("outerCore",      -H/2,          H/2,           CORE_NODES,       ['O'])
    # Upper axial regions (fuel assemblies)
    _add_zone("upperGasPlenum",  H/2,           H/2+UF,        UPPER_FGP_NODES,  ['I','O'])
    _add_zone("upperReflector",  H/2+UF,        H/2+UF+UR,     UPPER_REFL_NODES, ['I','O'])
    # Lower axial regions (fuel assemblies)
    _add_zone("lowerGasPlenum", -H/2-LF,       -H/2,           LOW_FGP_NODES,    ['I','O'])
    _add_zone("lowerReflector", -H/2-LF-LR,    -H/2-LF,       LOW_REFL_NODES,   ['I','O'])
    # Radial reflector (all axial levels)
    for zmin, zmax, nz in [
        (-H/2, H/2, CORE_NODES), (H/2, H/2+UF, UPPER_FGP_NODES),
        (H/2+UF, H/2+UF+UR, UPPER_REFL_NODES), (-H/2-LF, -H/2, LOW_FGP_NODES),
        (-H/2-LF-LR, -H/2-LF, LOW_REFL_NODES)
    ]:
        _add_zone("radialReflector", zmin, zmax, nz, ['R'])
    # Diagrid support structure (below all reflectors)
    _add_zone("diagrid",   -H/2-LF-LR-DG, -H/2-LF-LR, DIAGRID_NODES, ['I','O','R','C'])
    # Follower (control rod absorber segment below active core)
    for zmin, zmax, nz in [
        (-H/2-LF-LR, -H/2-LF, LOW_REFL_NODES),
        (-H/2-LF,    -H/2,    LOW_FGP_NODES),
        (-H/2,        H/2,    CORE_NODES),
    ]:
        _add_zone("follower", zmin, zmax, nz, ['C'])
    # Control rod absorber (above active core, rods inserted from the top)
    _add_zone("controlRod", H/2,     H/2+UF,    UPPER_FGP_NODES,  ['C'])
    _add_zone("controlRod", H/2+UF,  H/2+UF+UR, UPPER_REFL_NODES, ['C'])

    # Connect top/bottom faces of adjacent axial zones into internal interfaces
    nMesh.add_merge_patch_pairs(includeFacename=['Top','Bottom'], excludeFacename=['Wall'])
    # Weld side faces between neighbouring assemblies into internal interfaces
    nMesh.add_merge_patch_pairs(includeFacename=['Wall'], excludeFacename=['Top','Bottom'])

    ext_walls = nMesh.get_standalone_faces(
        includeFacename=['Wall'], excludeFacename=['Top','Bottom']
    )
    nMesh.merge_patches_with_name(
        name="wall",
        includeFacename=[f.name for f in ext_walls], patchType="wall"
    )
    nMesh.merge_patches_with_name(name="bottom", includeFacename=["diagridBottom_"])
    nMesh.merge_patches_with_name(
        name="top",
        includeFacename=["upperReflectorTop_", "controlRodTop_", "radialReflectorTop_"]
    )
