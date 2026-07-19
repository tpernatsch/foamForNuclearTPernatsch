import foamForNuclear as ffn
from foamForNuclear import mesh


def build_fluid_mesh():
    """
    Create a 2D segmented fluid channel with:
    - upstream (slip bottom)
    - middle (interface)
    - downstream (wall bottom)
    """

    fluid_mesh = mesh.BlockMesh(region="fluid")

    # --- Blocks ---
    fluidBlock2 = fluid_mesh.create_cube(
        "zone0", 0, 0, 0,
        1, 0.5, 0.4,
        nx=200, ny=41, nz=1,
        gradx=5, grady=16
    )

    fluidBlock3 = fluid_mesh.extrude_right(
        [fluidBlock2], "zone0",
        dx=2, nx=51
    )

    fluidBlock1 = fluid_mesh.extrude_left(
        [fluidBlock2], "zone0",
        dx=0.5, nx=81, gradx=0.2
    )

    # --- Boundaries ---

    # Inlet
    inlet = ffn.mesh.Face("inlet")
    inlet.add_sub_face(fluidBlock1.leftFace())

    # Outlet
    outlet = ffn.mesh.Face("outlet")
    outlet.add_sub_face(fluidBlock3.rightFace())

    # Top (merged)
    top = ffn.mesh.Face("top", boundaryType="wall")
    for blk in [fluidBlock1, fluidBlock2, fluidBlock3]:
        top.add_sub_face(blk.backFace())

    # Bottom downstream (wall)
    bottom = ffn.mesh.Face("bottom")
    bottom.add_sub_face(fluidBlock3.frontFace())

    # Bottom upstream (slip)
    slip_bottom = ffn.mesh.Face("slip-bottom")
    slip_bottom.add_sub_face(fluidBlock1.frontFace())

    # Interface (middle)
    interface = ffn.mesh.Face(
        "interface",
        boundaryType="mappedWall",
        inGroups=["wall"],
        extraParameters={
            "sampleRegion": "solid",
            "samplePatch": "top",
            "sampleMode": "nearestPatchFace"
        }
    )
    interface.add_sub_face(fluidBlock2.frontFace())

    # --- Add boundaries ---
    for b in [inlet, outlet, top, bottom, slip_bottom, interface]:
        fluid_mesh.add_boundary(b)

    return fluid_mesh


def build_solid_mesh():
    """
    Create a 2D solid block under the middle fluid section.
    """

    solid_mesh = mesh.BlockMesh(region="solid")

    solidBlock = solid_mesh.create_cube(
        "solid",
        0, -0.25, 0,
        1, 0, 0.4,
        nx=200, ny=41, nz=1,
        gradx=5, grady=0.0625
    )

    # --- Boundaries ---

    left = ffn.mesh.Face("left", boundaryType="wall")
    left.add_sub_face(solidBlock.leftFace())

    right = ffn.mesh.Face("right", boundaryType="wall")
    right.add_sub_face(solidBlock.rightFace())

    top = ffn.mesh.Face(
        "top",
        boundaryType="mappedWall",
        inGroups=["wall"],
        extraParameters={
            "sampleRegion": "fluid",
            "samplePatch": "interface",
            "sampleMode": "nearestPatchFace"
        }
    )
    top.add_sub_face(solidBlock.backFace())

    bottom = ffn.mesh.Face("bottom", boundaryType="wall")
    bottom.add_sub_face(solidBlock.frontFace())

    # --- Add boundaries ---
    for b in [left, right, top, bottom]:
        solid_mesh.add_boundary(b)

    return solid_mesh


def build_mesh():
    """
    Convenience wrapper returning both meshes.
    """
    fluid = build_fluid_mesh()
    solid = build_solid_mesh()
    return fluid, solid