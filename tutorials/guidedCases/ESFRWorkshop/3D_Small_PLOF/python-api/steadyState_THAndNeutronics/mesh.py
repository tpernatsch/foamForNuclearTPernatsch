"""
mesh.py — Computational meshes for the coupled TH + neutronics case
====================================================================

This case uses TWO separate meshes that occupy the same physical space
but serve different physics solvers:

  nMesh  — neutronics mesh (neutroRegion)
            Coarser hexagonal prism mesh; one volume per assembly per
            axial layer.  The neutron flux is homogenised over each volume.

  thMesh — thermal-hydraulics mesh (fluidRegion)
            Same geometry as nMesh but with internal baffles added between
            assemblies to enforce separate upward flow channels.

Fields are exchanged between the two meshes at every coupling iteration
(power density from neutronics to TH, temperatures from TH to neutronics).
This exchange is set up in caseSetup.py.

TWO MESH OPTIONS:  Option A (load pre-generated) or Option B (generate)
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
    nMesh  = mesh.PolyMesh(region="neutroRegion", srcpath="../meshes/polyMeshNeutro")
    thMesh = mesh.PolyMesh(region="fluidRegion",  srcpath="../meshes/polyMeshFluid")

else:
    # Shared geometry parameters
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

    H  = CORE_HEIGHT; UF = UPPER_FGP_HEIGHT; UR = UPPER_REFL_HEIGHT
    LF = LOW_FGP_HEIGHT; LR = LOW_REFL_HEIGHT; DG = DIAGRID_HEIGHT

    def _build_mesh(region_name, is_fluid):
        m = mesh.BlockMesh(region=region_name)

        def _add(name, zmin, zmax, nz, elements):
            m.lattice_placement(
                funcElementGenerator=lambda n, x, y: m.create_hexagon_prism_along_z(
                    name=name, zmin=zmin, zmax=zmax,
                    pitch=flat_to_flat, x=x, y=y,
                    nr=nr, nt=nt, nz=nz, isAddAllBC=True
                ),
                lattice=LATTICE, latticeType="hexagon",
                nx=lattice_n, ny=lattice_n, pitch=pitch, elementsToPlace=elements,
            )

        _add("innerCore",      -H/2, H/2, CORE_NODES, ['I'])
        _add("outerCore",      -H/2, H/2, CORE_NODES, ['O'])
        _add("upperGasPlenum",  H/2,  H/2+UF,      UPPER_FGP_NODES,  ['I','O'])
        _add("upperReflector",  H/2+UF, H/2+UF+UR, UPPER_REFL_NODES, ['I','O'])
        _add("lowerGasPlenum", -H/2-LF,    -H/2,   LOW_FGP_NODES,    ['I','O'])
        _add("lowerReflector", -H/2-LF-LR, -H/2-LF, LOW_REFL_NODES,  ['I','O'])
        for zmin, zmax, nz in [
            (-H/2, H/2, CORE_NODES), (H/2, H/2+UF, UPPER_FGP_NODES),
            (H/2+UF, H/2+UF+UR, UPPER_REFL_NODES), (-H/2-LF, -H/2, LOW_FGP_NODES),
            (-H/2-LF-LR, -H/2-LF, LOW_REFL_NODES)
        ]:
            _add("radialReflector", zmin, zmax, nz, ['R'])
        _add("diagrid", -H/2-LF-LR-DG, -H/2-LF-LR, DIAGRID_NODES, ['I','O','R','C'])
        for zmin, zmax, nz in [
            (-H/2-LF-LR, -H/2-LF, LOW_REFL_NODES),
            (-H/2-LF,    -H/2,    LOW_FGP_NODES),
            (-H/2,        H/2,    CORE_NODES),
        ]:
            _add("follower", zmin, zmax, nz, ['C'])
        _add("controlRod", H/2,    H/2+UF,    UPPER_FGP_NODES,  ['C'])
        _add("controlRod", H/2+UF, H/2+UF+UR, UPPER_REFL_NODES, ['C'])

        m.add_merge_patch_pairs(includeFacename=['Top','Bottom'], excludeFacename=['Wall'])
        if is_fluid:
            m.add_baffles(includeFacename=['Wall'], excludeFacename=['Top','Bottom'])
        else:
            m.add_merge_patch_pairs(includeFacename=['Wall'], excludeFacename=['Top','Bottom'])

        ext_walls = m.get_standalone_faces(includeFacename=['Wall'], excludeFacename=['Top','Bottom'])
        m.merge_patches_with_name(name="wall", includeFacename=[f.name for f in ext_walls], patchType="wall")
        m.merge_patches_with_name(name="bottom", includeFacename=["diagridBottom_"])
        m.merge_patches_with_name(name="top", includeFacename=["upperReflectorTop_", "controlRodTop_", "radialReflectorTop_"])
        return m

    nMesh  = _build_mesh("neutroRegion", is_fluid=False)
    thMesh = _build_mesh("fluidRegion",  is_fluid=True)
