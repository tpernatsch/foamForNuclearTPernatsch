from __future__ import annotations

import keyword
from typing import Literal, Optional
import foamForNuclear.mesh as mesh

import attrs

@attrs.define(init=False, slots=False, kw_only=True)
class Rod1DBlockMesh(mesh.BlockMesh):
    scale: float = 1.0
    wedge_angle: float = 0.5
    fuel_length: float
    layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad"
    fuel_ri: float
    fuel_ro: float
    fuel_nr: float
    fuel_gradng_nr: float = 1
    clad_ri: float
    clad_ro: float
    clad_nr: int
    plenum_length: float = 0.0
    slices: int = 1
    fuel_name: str = "fuel"
    clad_name: str = "clad"
    offset: float = 0.0

    def __init__(
        self,
        *,
        scale: float = 1.0,
        wedge_angle: float = 0.5,
        fuel_length: float,
        fuel_ri: Optional[float] = None,
        fuel_ro: Optional[float] = None,
        clad_ri: Optional[float] = None,
        clad_ro: Optional[float] = None,
        plenum_length: float = 0.0,
        fuel_name: str = "fuel",
        clad_name: str = "cladding",
        offset: float = 0.0,
        fuel_nr: Optional[int] = None,
        fuel_grading_nr: Optional[float] = 1,
        clad_nr: Optional[int] = None,
        slices: int = 1,
        layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad",
    ):
        # init BlockMesh first
        super().__init__(scale=scale)

        # then set attrs fields
        self.layout = layout
        self.fuel_length = fuel_length
        self.fuel_ri = fuel_ri
        self.fuel_ro = fuel_ro
        self.fuel_nr = fuel_nr
        self.fuel_grading_nr = fuel_grading_nr
        self.clad_ri = clad_ri
        self.clad_ro = clad_ro
        self.clad_nr = clad_nr
        self.plenum_length = plenum_length
        self.slices = slices
        self.fuel_name = fuel_name
        self.clad_name = clad_name
        self.offset = offset
        self.wedge_angle = wedge_angle
        self.scale = scale

def rod_1d(
    *,
    scale: float = 1.0,
    wedge_angle: float = 0.5,
    fuel_length: float,
    fuel_ri: Optional[float] = None,
    fuel_ro: Optional[float] = None,
    clad_ri: Optional[float] = None,
    clad_ro: Optional[float] = None,
    plenum_length: float = 0.0,
    fuel_name: str = "fuel",
    clad_name: str = "cladding",
    offset: float = 0.0,
    fuel_nr: Optional[int] = None,
    fuel_grading_nr: Optional[float] = 1,
    clad_nr: Optional[int] = None,
    slices: int = 1,
    neighbour_region: str = "region0",
    layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad",
    updateAMI: bool = False
) -> Rod1DBlockMesh:
    """
    Build a 1D wedge rod mesh with fuel and/or cladding and an optional top plenum
    extrusion.

    Parameters
    ----------
    scale : float, optional
        Global mesh scale factor (passed to BlockMesh).
    wedge_angle : float, optional
        Wedge angle (radians or whatever the underlying mesh expects).
    fuel_length : float
        Axial length of the fuel / cladding block (before plenum).
    fuel_ri, fuel_ro : float, optional
        Fuel inner/outer radii. Required if layout uses fuel
        ("fuel_only", "fuel_and_clad").
    clad_ri, clad_ro : float, optional
        Cladding inner/outer radii. Required if layout uses cladding
        ("clad_only", "fuel_and_clad").
    plenum_length : float, optional
        Extra axial length to extrude at the top:
        - for layout="fuel_and_clad": top of cladding only (OFFBEAT behaviour),
        - for layout="fuel_only": top of fuel,
        - for layout="clad_only": top of cladding.
        (0 ⇒ no plenum.)
    fuel_name : str, optional
        Region/name for the fuel wedge.
    clad_name : str, optional
        Region/name for the cladding wedge.
    offset : float, optional
        Axial offset for the wedge origin (passed through to create_wedge).
    fuel_nr, clad_nr : int, optional
        Radial cell counts for fuel/clad wedges, required depending on layout.
    fuel_grading_nr: Optionl[float] = 1
        Radial cell grading: a number smaller than 1 makes the fuel mesh progressively finer towards the rim.        
    slices : int, optional
        Axial cell count (z-direction) for the base blocks (before plenum).
    neighbour_region : str, optional
        Value written in the region-coupled BCs (neighbourRegion key) when
        layout="fuel_and_clad".
    layout : {"fuel_and_clad", "fuel_only", "clad_only"}, optional
        Controls which regions are present:
        - "fuel_and_clad" (default): fuel + cladding, with a regionCoupled
          interface between fuelOuter and cladInner (OFFBEAT-style).
        - "fuel_only": only the fuel region.
        - "clad_only": only the cladding region.

    Returns
    -------
    mesh.BlockMesh
        The assembled mesh object with boundaries added.
    """

    has_fuel = layout in ("fuel_only", "fuel_and_clad")
    has_clad = layout in ("clad_only", "fuel_and_clad")

    # ---- basic global checks ----
    if fuel_length <= 0.0:
        raise ValueError("fuel_length must be > 0.")
    if plenum_length < 0.0:
        raise ValueError("plenum_length cannot be negative.")
    if slices <= 0:
        raise ValueError("slices must be a positive integer.")

    # ---- layout-dependent input checks ----
    if has_fuel:
        if fuel_ri is None or fuel_ro is None or fuel_nr is None:
            raise ValueError(
                "layout includes fuel: fuel_ri, fuel_ro and fuel_nr must be provided."
            )
        if not (fuel_ri >= 0.0 and fuel_ro > fuel_ri):
            raise ValueError("Fuel radii must satisfy: fuel_ri >= 0 and fuel_ro > fuel_ri.")
        if fuel_nr <= 0:
            raise ValueError("fuel_nr must be a positive integer.")

    if has_clad:
        if clad_ri is None or clad_ro is None or clad_nr is None:
            raise ValueError(
                "layout includes cladding: clad_ri, clad_ro and clad_nr must be provided."
            )
        if not (clad_ri >= 0.0 and clad_ro > clad_ri):
            raise ValueError("Cladding radii must satisfy: clad_ri >= 0 and clad_ro > clad_ri.")
        if clad_nr <= 0:
            raise ValueError("clad_nr must be a positive integer.")

    # extra consistency check for fuel+clad
    if has_fuel and has_clad:
        if not (clad_ri >= fuel_ro):
            raise ValueError(
                "Expected clad_ri >= fuel_ro for fuel+clad layout "
                "(gap >= 0, no overlap)."
            )

    # ---- mesh + primary blocks ----
    bm = Rod1DBlockMesh(
        scale=scale,
        wedge_angle=wedge_angle,
        fuel_length=fuel_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        fuel_nr=fuel_nr,
        fuel_grading_nr=fuel_grading_nr,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        clad_nr=clad_nr,
        plenum_length=plenum_length,
        fuel_name=fuel_name,
        clad_name=clad_name,
        offset=offset,
        slices=slices,
        layout=layout,
    )

    fuel_block = None
    clad_block = None
    fuel_block_top = None
    clad_block_top = None

    if has_fuel:
        fuel_block = bm.create_wedge(
            fuel_name,
            fuel_ri,
            fuel_ro,
            offset,
            fuel_length,
            wedgeAngle=wedge_angle,
            nr=fuel_nr,
            nz=slices,
            gradr=fuel_grading_nr,
        )

    if has_clad:
        clad_block = bm.create_wedge(
            clad_name,
            clad_ri,
            clad_ro,
            offset,
            fuel_length,
            wedgeAngle=wedge_angle,
            nr=clad_nr,
            nz=slices,
        )

    # ---- plenum extrusion ----
    if plenum_length > 0.0:
        if has_clad and layout == "fuel_and_clad":
            # reproduce your original behaviour: plenum on cladding only
            clad_block_top = bm.extrude_top([clad_block], clad_name, dz=plenum_length, nz=1)
        elif has_fuel and not has_clad:
            # fuel-only case: extrude the fuel block
            fuel_block_top = bm.extrude_top([fuel_block], fuel_name, dz=plenum_length, nz=1)
        elif has_clad and not has_fuel:
            # clad-only case: extrude the cladding block
            clad_block_top = bm.extrude_top([clad_block], clad_name, dz=plenum_length, nz=1)

    # ---- wedge faces (front/back) ----
    if has_fuel:
        bc_fuel_front = mesh.Face(name="fuelFront", boundaryType="wedge")
        bc_fuel_front.add_sub_face(fuel_block.frontFace())
        if fuel_block_top is not None:
            bc_fuel_front.add_sub_face(fuel_block_top.frontFace())
        bm.add_boundary(bc_fuel_front)

        bc_fuel_back = mesh.Face(name="fuelBack", boundaryType="wedge")
        bc_fuel_back.add_sub_face(fuel_block.backFace())
        if fuel_block_top is not None:
            bc_fuel_back.add_sub_face(fuel_block_top.backFace())
        bm.add_boundary(bc_fuel_back)

    if has_clad:
        bc_clad_front = mesh.Face(name="cladFront", boundaryType="wedge")
        bc_clad_front.add_sub_face(clad_block.frontFace())
        if clad_block_top is not None:
            bc_clad_front.add_sub_face(clad_block_top.frontFace())
        bm.add_boundary(bc_clad_front)

        bc_clad_back = mesh.Face(name="cladBack", boundaryType="wedge")
        bc_clad_back.add_sub_face(clad_block.backFace())
        if clad_block_top is not None:
            bc_clad_back.add_sub_face(clad_block_top.backFace())
        bm.add_boundary(bc_clad_back)

    # ---- radial boundaries ----
    update_ami_str = "true" if updateAMI else "false"
    if has_fuel and has_clad:
        # region-coupled interface: fuelOuter <-> cladInner
        bc_fuel_outer = mesh.Face(
            name="fuelOuter",
            boundaryType="regionCoupledOFFBEAT",
            extraParameters={
                "neighbourPatch": "cladInner",
                "neighbourRegion": neighbour_region,
                "owner": "true",
                "updateAMI": f"{update_ami_str}",
            },
        )
        bc_fuel_outer.add_sub_face(fuel_block.rightFace())
        bm.add_boundary(bc_fuel_outer)

        bc_clad_inner = mesh.Face(
            name="cladInner",
            boundaryType="regionCoupledOFFBEAT",
            extraParameters={
                "neighbourPatch": "fuelOuter",
                "neighbourRegion": neighbour_region,
                "owner": "false",
                "updateAMI": f"{update_ami_str}",
            },
        )
        bc_clad_inner.add_sub_face(clad_block.leftFace())
        if clad_block_top is not None:
            bc_clad_inner.add_sub_face(clad_block_top.leftFace())
        bm.add_boundary(bc_clad_inner)

        # external cladding outer wall
        bc_clad_outer = mesh.Face(name="cladOuter")
        bc_clad_outer.add_sub_face(clad_block.rightFace())
        if clad_block_top is not None:
            bc_clad_outer.add_sub_face(clad_block_top.rightFace())
        bm.add_boundary(bc_clad_outer)

    elif has_fuel and not has_clad:
        # only fuel: outer wall is a "normal" patch
        bc_fuel_outer = mesh.Face(name="fuelOuter")
        bc_fuel_outer.add_sub_face(fuel_block.rightFace())
        if fuel_block_top is not None:
            bc_fuel_outer.add_sub_face(fuel_block_top.rightFace())
        bm.add_boundary(bc_fuel_outer)

    elif has_clad and not has_fuel:
        # clad-only: inner and outer are both "real" boundaries
        bc_clad_inner = mesh.Face(name="cladInner")
        bc_clad_inner.add_sub_face(clad_block.leftFace())
        if clad_block_top is not None:
            bc_clad_inner.add_sub_face(clad_block_top.leftFace())
        bm.add_boundary(bc_clad_inner)

        bc_clad_outer = mesh.Face(name="cladOuter")
        bc_clad_outer.add_sub_face(clad_block.rightFace())
        if clad_block_top is not None:
            bc_clad_outer.add_sub_face(clad_block_top.rightFace())
        bm.add_boundary(bc_clad_outer)

    # ---- caps (top/bottom) ----
    cap_faces: list[tuple[str, mesh.Face]] = []

    if has_fuel:
        cap_faces.append(("fuelBottom", fuel_block.bottomFace()))
        top_face = fuel_block_top.topFace() if fuel_block_top is not None else fuel_block.topFace()
        cap_faces.append(("fuelTop", top_face))

    if has_clad:
        cap_faces.append(("cladBottom", clad_block.bottomFace()))
        top_face = clad_block_top.topFace() if clad_block_top is not None else clad_block.topFace()
        cap_faces.append(("cladTop", top_face))

    for bc_name, subface in cap_faces:
        bc_empty = mesh.Face(name=bc_name, boundaryType="empty")
        bc_empty.add_sub_face(subface)
        bm.add_boundary(bc_empty)

    return bm