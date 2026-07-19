from __future__ import annotations

from typing import Optional
import attrs
import foamForNuclear.mesh as mesh


@attrs.define(init=False, slots=False, kw_only=True)
class Sphere1DBlockMesh(mesh.BlockMesh):
    """
    BlockMesh for a 1D spherical TRISO-style geometry with concentric layers.

    Attributes
    ----------
    wedge_angle : float
        Wedge opening angle in degrees.
    inner_radius : float
        Inner radius of the innermost modeled layer.
    layers : list of (str, float, int)
        Layer definitions: (zone_name, outer_radius, nr_radial).
    """
    wedge_angle: float = 0.25
    inner_radius: float = 0.0
    layers: list

    def __init__(
        self,
        *,
        wedge_angle: float = 0.25,
        scale: float = 1.0,
        inner_radius: float = 0.0,
        layers: list,
    ):
        super().__init__(scale=scale)
        self.wedge_angle = wedge_angle
        self.inner_radius = inner_radius
        self.layers = layers


def sphere_1d(
    *,
    wedge_angle: float = 0.25,
    scale: float = 1.0,
    inner_radius: float = 0.0,
    layers: list[tuple[str, float, int]],
    sample_region: str | None = None,
) -> Sphere1DBlockMesh:
    """
    Build a 1D spherical (TRISO-style) mesh with concentric layers.

    Parameters
    ----------
    wedge_angle : float
        Wedge opening angle in degrees (default 0.25).
    scale : float
        Global scale factor for the mesh (e.g. 1e-6 to input radii in μm).
    inner_radius : float
        Inner radius of the first layer. If > 0, an ``inner`` patch is created.
        If 0, no inner patch (solid sphere starting at origin).
    layers : list of (name, outer_radius, nr_radial)
        Ordered innermost → outermost. Each tuple:

        - ``name``: cell zone / material name (e.g. ``"SiC"``, ``"IPyC"``).
        - ``outer_radius``: outer boundary of this layer (in scale units).
        - ``nr_radial``: number of radial cells in this layer.
    sample_region : str or None
        When set, the outer patch is created as ``mappedPatch`` sampling from
        ``sample_region`` in ``nearestCell`` mode.  Required for
        ``mappedFixedValue`` BCs (e.g. to read T from a macro compact mesh).
        When ``None`` (default), the outer patch is a plain ``patch``.

    Returns
    -------
    Sphere1DBlockMesh
        BlockMesh with ``front``, ``back``, ``top``, ``bottom`` (wedge) patches,
        an ``inner`` patch (only when ``inner_radius > 0``), and an ``outer``
        patch at the outermost surface.
    """
    if not layers:
        raise ValueError("At least one layer must be provided.")

    bm = Sphere1DBlockMesh(
        wedge_angle=wedge_angle,
        scale=scale,
        inner_radius=inner_radius,
        layers=layers,
    )
    bm.isMergeCoincidentPoints = True

    blocks = []
    r_in = inner_radius
    for (name, r_out, nr) in layers:
        block = bm.create_sphere_1D(
            name=name,
            innerRadius=r_in,
            outerRadius=r_out,
            opening=wedge_angle,
            nr=nr,
        )
        blocks.append(block)
        r_in = r_out

    # 4 wedge boundary patches (matching sphereMaker.py patch names)
    for patch_name, method in [
        ("front",  "frontFace"),
        ("back",   "backFace"),
        ("top",    "topFace"),
        ("bottom", "bottomFace"),
    ]:
        bc = mesh.Face(name=patch_name, boundaryType="wedge")
        for block in blocks:
            bc.add_sub_face(getattr(block, method)())
        bm.add_boundary(bc)

    # Inner patch (only when the inner surface is exposed)
    if inner_radius > 0.0:
        bc_inner = mesh.Face(name="inner", boundaryType="patch")
        bc_inner.add_sub_face(blocks[0].leftFace())
        bm.add_boundary(bc_inner)

    # Outer patch (outermost surface)
    if sample_region is not None:
        bc_outer = mesh.Face(
            name="outer",
            boundaryType="mappedPatch",
            inGroups=["mappedPatch"],
            extraParameters={
                "sampleMode":   "nearestCell",
                "sampleRegion": sample_region,
                "offsetMode":   "uniform",
                "offset":       "(0 0 0)",
            },
        )
    else:
        bc_outer = mesh.Face(name="outer", boundaryType="patch")
    bc_outer.add_sub_face(blocks[-1].rightFace())
    bm.add_boundary(bc_outer)

    return bm
