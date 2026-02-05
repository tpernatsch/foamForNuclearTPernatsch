from __future__ import annotations

from typing import Literal, Optional

import attrs
import foamForNuclear.mesh as mesh


@attrs.define(init=False, slots=False, kw_only=True)
class Rod2DRZBlockMesh(mesh.BlockMesh):
    scale: float = 1.0
    wedge_angle: float = 0.5

    fuel_length: float
    plenum_length: float = 0.0

    layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad"

    fuel_ri: float
    fuel_ro: float
    fuel_nr: int
    fuel_grading_nr: float = 1.0
    fuel_nz: int = 1

    clad_ri: float
    clad_ro: float
    clad_nr: int
    clad_nz: int = 1  # for the *whole* clad block (fuel_length + plenum_length)

    fuel_name: str = "fuel"
    clad_name: str = "cladding"
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
        fuel_grading_nr: float = 1.0,
        clad_nr: Optional[int] = None,
        fuel_nz: int = 1,
        clad_nz: int = 1,
        layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad",
    ):
        super().__init__(scale=scale)

        self.scale = scale
        self.wedge_angle = wedge_angle

        self.fuel_length = fuel_length
        self.plenum_length = plenum_length
        self.layout = layout

        self.fuel_ri = fuel_ri
        self.fuel_ro = fuel_ro
        self.fuel_nr = fuel_nr
        self.fuel_grading_nr = fuel_grading_nr
        self.fuel_nz = fuel_nz

        self.clad_ri = clad_ri
        self.clad_ro = clad_ro
        self.clad_nr = clad_nr
        self.clad_nz = clad_nz

        self.fuel_name = fuel_name
        self.clad_name = clad_name
        self.offset = offset


def rod_2d_rz(
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
    fuel_grading_nr: float = 1.0,
    clad_nr: Optional[int] = None,
    fuel_nz: int = 1,
    clad_nz: int = 1,
    neighbour_region: str = "region0",
    layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad",
    updateAMI: bool = False
) -> Rod2DRZBlockMesh:
    """
    Build a 2D r-z wedge rod mesh with fuel and/or cladding.

    Differences vs rod_1d:
      - no `slices`; you provide `fuel_nz` and `clad_nz`
      - cladding can be a single block of height (fuel_length + plenum_length)
        divided into `clad_nz`
      - top and bottom faces are patches (not `empty`)
    """

    has_fuel = layout in ("fuel_only", "fuel_and_clad")
    has_clad = layout in ("clad_only", "fuel_and_clad")

    if fuel_length <= 0.0:
        raise ValueError("fuel_length must be > 0.")
    if plenum_length < 0.0:
        raise ValueError("plenum_length cannot be negative.")
    if fuel_nz <= 0:
        raise ValueError("fuel_nz must be a positive integer.")
    if clad_nz <= 0:
        raise ValueError("clad_nz must be a positive integer.")

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

    if has_fuel and has_clad:
        if not (clad_ri >= fuel_ro):
            raise ValueError(
                "Expected clad_ri >= fuel_ro for fuel+clad layout (gap >= 0, no overlap)."
            )

    bm = Rod2DRZBlockMesh(
        scale=scale,
        wedge_angle=wedge_angle,
        fuel_length=fuel_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        fuel_nr=fuel_nr,
        fuel_grading_nr=fuel_grading_nr,
        fuel_nz=fuel_nz,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        clad_nr=clad_nr,
        clad_nz=clad_nz,
        plenum_length=plenum_length,
        fuel_name=fuel_name,
        clad_name=clad_name,
        offset=offset,
        layout=layout,
    )

    fuel_block = None
    clad_block = None

    if has_fuel:
        fuel_block = bm.create_wedge(
            fuel_name,
            fuel_ri,
            fuel_ro,
            offset,
            fuel_length,
            wedgeAngle=wedge_angle,
            nr=fuel_nr,
            nz=fuel_nz,
            gradr=fuel_grading_nr,
        )

    if has_clad:
        clad_height = fuel_length + plenum_length
        clad_block = bm.create_wedge(
            clad_name,
            clad_ri,
            clad_ro,
            offset,
            clad_height,
            wedgeAngle=wedge_angle,
            nr=clad_nr,
            nz=clad_nz,
        )

    # ---- wedge faces (front/back) ----
    if has_fuel:
        bc_fuel_front = mesh.Face(name="fuelFront", boundaryType="wedge")
        bc_fuel_front.add_sub_face(fuel_block.frontFace())
        bm.add_boundary(bc_fuel_front)

        bc_fuel_back = mesh.Face(name="fuelBack", boundaryType="wedge")
        bc_fuel_back.add_sub_face(fuel_block.backFace())
        bm.add_boundary(bc_fuel_back)

    if has_clad:
        bc_clad_front = mesh.Face(name="cladFront", boundaryType="wedge")
        bc_clad_front.add_sub_face(clad_block.frontFace())
        bm.add_boundary(bc_clad_front)

        bc_clad_back = mesh.Face(name="cladBack", boundaryType="wedge")
        bc_clad_back.add_sub_face(clad_block.backFace())
        bm.add_boundary(bc_clad_back)

    # ---- radial boundaries ----    
    update_ami_str = "true" if updateAMI else "false"
    if has_fuel and has_clad:
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
        bm.add_boundary(bc_clad_inner)

        bc_clad_outer = mesh.Face(name="cladOuter")
        bc_clad_outer.add_sub_face(clad_block.rightFace())
        bm.add_boundary(bc_clad_outer)

    elif has_fuel and not has_clad:
        bc_fuel_outer = mesh.Face(name="fuelOuter")
        bc_fuel_outer.add_sub_face(fuel_block.rightFace())
        bm.add_boundary(bc_fuel_outer)

    elif has_clad and not has_fuel:
        bc_clad_inner = mesh.Face(name="cladInner")
        bc_clad_inner.add_sub_face(clad_block.leftFace())
        bm.add_boundary(bc_clad_inner)

        bc_clad_outer = mesh.Face(name="cladOuter")
        bc_clad_outer.add_sub_face(clad_block.rightFace())
        bm.add_boundary(bc_clad_outer)

    # ---- caps (top/bottom) as PATCH (not empty) ----
    cap_faces: list[tuple[str, mesh.Face]] = []

    if has_fuel:
        cap_faces.append(("fuelBottom", fuel_block.bottomFace()))
        cap_faces.append(("fuelTop", fuel_block.topFace()))

    if has_clad:
        cap_faces.append(("cladBottom", clad_block.bottomFace()))
        cap_faces.append(("cladTop", clad_block.topFace()))

    for bc_name, subface in cap_faces:
        bc_patch = mesh.Face(name=bc_name, boundaryType="patch")
        bc_patch.add_sub_face(subface)
        bm.add_boundary(bc_patch)

    return bm