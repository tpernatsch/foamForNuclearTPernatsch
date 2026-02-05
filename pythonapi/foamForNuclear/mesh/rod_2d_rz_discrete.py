from __future__ import annotations

from typing import Literal, Optional
import attrs
import foamForNuclear.mesh as mesh

import math

@attrs.define(init=False, slots=False, kw_only=True)
class Rod2DRZDiscreteBlockMesh(mesh.BlockMesh):
    scale: float = 1.0
    wedge_angle: float = 0.5

    # layout
    layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad"

    # radii
    fuel_ri: float = 0.0
    fuel_ro: float
    clad_ri: float
    clad_ro: float

    # discretization
    fuel_nr: int
    fuel_grading_nr: float = 1.0
    pellet_count: int
    pellet_nz: int
    clad_nr: int
    clad_nz: int

    # stack
    pellet_height: float
    plenum_length: float = 0.0
    offset: float = 0.0

    # pellet geometry
    pellet_type: Literal["flat", "dished", "chamfered", "dished_and_chamfered"] = "flat"

    # dishing
    dish_radius_curvature: Optional[float] = None
    dish_outer_radius: Optional[float] = None  # radius that separates dished inner ring from outer ring

    # chamfer
    chamfer_height: Optional[float] = None
    chamfer_start_radius: Optional[float] = None  # where chamfer ring begins (inner radius of chamfered ring)

    # names
    fuel_name: str = "fuel"
    clad_name: str = "cladding"

    def __init__(
        self,
        *,
        scale: float = 1.0,
        wedge_angle: float = 0.5,
        layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad",

        fuel_ri: float = 0.0,
        fuel_ro: float,
        clad_ri: Optional[float] = None,
        clad_ro: Optional[float] = None,

        fuel_nr: int,
        fuel_grading_nr: float = 1.0,
        pellet_count: int,
        pellet_nz: int,
        clad_nr: Optional[int] = None,
        clad_nz: Optional[int] = None,

        pellet_height: float,
        plenum_length: float = 0.0,
        offset: float = 0.0,

        pellet_type: Literal["flat", "dished", "chamfered", "dished_and_chamfered"] = "flat",
        dish_radius_curvature: Optional[float] = None,
        dish_outer_radius: Optional[float] = None,
        chamfer_height: Optional[float] = None,
        chamfer_start_radius: Optional[float] = None,

        fuel_name: str = "fuel",
        clad_name: str = "cladding",
    ):
        super().__init__(scale=scale)

        self.scale = scale
        self.wedge_angle = wedge_angle
        self.layout = layout

        self.fuel_ri = fuel_ri
        self.fuel_ro = fuel_ro
        self.clad_ri = clad_ri
        self.clad_ro = clad_ro

        self.fuel_nr = fuel_nr
        self.fuel_grading_nr = fuel_grading_nr
        self.pellet_count = pellet_count
        self.pellet_nz = pellet_nz
        self.clad_nr = clad_nr
        self.clad_nz = clad_nz

        self.pellet_height = pellet_height
        self.plenum_length = plenum_length
        self.offset = offset

        self.pellet_type = pellet_type
        self.dish_radius_curvature = dish_radius_curvature
        self.dish_outer_radius = dish_outer_radius
        self.chamfer_height = chamfer_height
        self.chamfer_start_radius = chamfer_start_radius

        self.fuel_name = fuel_name
        self.clad_name = clad_name


def rod_2d_rz_discrete(
    *,
    scale: float = 1.0,
    wedge_angle: float = 0.5,

    # geometry / radii
    fuel_ri: float = 0.0,
    fuel_ro: float,
    clad_ri: Optional[float] = None,
    clad_ro: Optional[float] = None,

    # discretization
    fuel_nr: int,
    fuel_grading_nr: float = 1.0,
    pellet_count: int,
    pellet_nz: int,
    clad_nr: Optional[int] = None,
    clad_nz: Optional[int] = None,

    # stack
    pellet_height: float,
    plenum_length: float = 0.0,
    offset: float = 0.0,

    # pellet shape
    pellet_type: Literal["flat", "dished", "chamfered", "dished_and_chamfered"] = "flat",

    # dishing (spherical cap)
    dish_radius_curvature: Optional[float] = None,   # R
    dish_outer_radius: Optional[float] = None,       # a

    # chamfer
    chamfer_height: Optional[float] = None,
    chamfer_start_radius: Optional[float] = None,    # start of chamfer ring (cald also "land" outer radius)

    # coupling / names
    neighbour_region: str = "region0",
    fuel_name: str = "fuel",
    clad_name: str = "cladding",
    layout: Literal["fuel_and_clad", "fuel_only", "clad_only"] = "fuel_and_clad",
    updateAMI: bool = False,

    # internal merging behaviour
    merge_radial_rings: bool = True,                 # merge ring interfaces inside each pellet
) -> "Rod2DRZDiscreteBlockMesh":
    """
    2D r-z axisymmetric rod with *discrete pellets* (stack of wedges).

    Strategy:
      - Create rings per pellet according to pellet_type.
      - Create explicit temporary patches for all internal interfaces.
      - Add deterministic mergePatchPairs by name (cheap, no overlap search).
      - Optionally aggregate/cleanup boundaries.

    Requires in BlockMesh:
      - create_wedge(...)
      - create_wedge_chamfered(...)
      - create_wedge_dished(..., dishRadiusCurvature=R, dishOuterRadius=a, ...)
      - mergePatchPairsByName(facename1, facename2)
      - mergeBoundaryFaces(...)
    """
    has_fuel = layout in ("fuel_only", "fuel_and_clad")
    has_clad = layout in ("clad_only", "fuel_and_clad")

    # -----------------------------
    # helpers
    # -----------------------------

    def _nr_split(r_lo: float, r_hi: float) -> int:
        """Stable radial cell split by thickness fraction (never returns 0)."""
        denom = (fuel_ro - fuel_ri)
        if denom <= 0:
            return max(1, fuel_nr)
        frac = (r_hi - r_lo) / denom
        return max(1, int(round(fuel_nr * frac)))

    def _dished_recess_at_inner_radius(R: float, a: float, ri: float) -> float:
        # dz_inner = sqrt(R^2 - ri^2) - sqrt(R^2 - a^2)
        return math.sqrt(R * R - ri * ri) - math.sqrt(R * R - a * a)

    # -----------------------------
    # basic checks
    # -----------------------------
    if pellet_count <= 0:
        raise ValueError("pellet_count must be >= 1.")
    if pellet_nz <= 0:
        raise ValueError("pellet_nz must be >= 1.")
    if pellet_height <= 0.0:
        raise ValueError("pellet_height must be > 0.")
    if plenum_length < 0.0:
        raise ValueError("plenum_length must be >= 0.")
    if fuel_nr <= 0:
        raise ValueError("fuel_nr must be >= 1.")
    if fuel_grading_nr <= 0:
        raise ValueError("fuel_grading_nr must be > 0.")

    if has_fuel:
        if fuel_ri < 0.0:
            raise ValueError("fuel_ri must be >= 0.")
        if fuel_ro <= fuel_ri:
            raise ValueError("Fuel radii must satisfy fuel_ro > fuel_ri.")
    if has_clad:
        if clad_ri is None or clad_ro is None or clad_nr is None or clad_nz is None:
            raise ValueError(
                "layout includes cladding: clad_ri, clad_ro, clad_nr, clad_nz are required."
            )
        if clad_ri < 0.0 or clad_ro <= clad_ri:
            raise ValueError("Clad radii must satisfy clad_ro > clad_ri >= 0.")
        if clad_nr <= 0 or clad_nz <= 0:
            raise ValueError("clad_nr and clad_nz must be >= 1.")
    if has_fuel and has_clad and clad_ri < fuel_ro:
        raise ValueError("Expected clad_ri >= fuel_ro for fuel+clad (gap >= 0, no overlap).")

    # -----------------------------
    # pellet-type checks (curvature-aware)
    # -----------------------------
    if pellet_type == "flat":
        pass

    elif pellet_type == "dished":
        if dish_radius_curvature is None or dish_outer_radius is None:
            raise ValueError("pellet_type='dished' requires dish_radius_curvature and dish_outer_radius.")
        R = float(dish_radius_curvature)
        a = float(dish_outer_radius)
        ri = float(fuel_ri)
        if R <= 0.0:
            raise ValueError("dish_radius_curvature (R) must be > 0.")
        if not (ri <= a < fuel_ro):
            raise ValueError("dish_outer_radius must satisfy fuel_ri <= dish_outer_radius < fuel_ro.")
        if R <= a:
            raise ValueError("Need dish_radius_curvature (R) > dish_outer_radius (a).")
        dz_inner = _dished_recess_at_inner_radius(R, a, ri)
        if dz_inner <= 0.0:
            raise ValueError("Dish geometry yields zero/negative recess at inner radius (check R, a, fuel_ri).")
        if dz_inner >= 0.5 * pellet_height:
            raise ValueError(
                f"Dish too deep for pellet height: dz_inner={dz_inner:g} "
                f"must be < pellet_height/2={0.5*pellet_height:g}."
            )

    elif pellet_type == "chamfered":
        if chamfer_height is None or chamfer_start_radius is None:
            raise ValueError("pellet_type='chamfered' requires chamfer_height and chamfer_start_radius.")
        if chamfer_height <= 0.0 or chamfer_height >= 0.5 * pellet_height:
            raise ValueError("chamfer_height must satisfy 0 < chamfer_height < pellet_height/2.")
        if not (fuel_ri <= chamfer_start_radius < fuel_ro):
            raise ValueError("chamfer_start_radius must satisfy fuel_ri <= chamfer_start_radius < fuel_ro.")

    elif pellet_type == "dished_and_chamfered":
        if dish_radius_curvature is None or dish_outer_radius is None or chamfer_start_radius is None or chamfer_height is None:
            raise ValueError(
                "pellet_type='dished_and_chamfered' requires dish_radius_curvature, dish_outer_radius, "
                "chamfer_start_radius, chamfer_height."
            )
        R = float(dish_radius_curvature)
        a = float(dish_outer_radius)
        ri = float(fuel_ri)
        if R <= 0.0:
            raise ValueError("dish_radius_curvature (R) must be > 0.")
        if not (ri <= a < chamfer_start_radius < fuel_ro):
            raise ValueError("Expected fuel_ri <= dish_outer_radius < chamfer_start_radius < fuel_ro.")
        if R <= a:
            raise ValueError("Need dish_radius_curvature (R) > dish_outer_radius (a).")
        dz_inner = _dished_recess_at_inner_radius(R, a, ri)
        if dz_inner <= 0.0:
            raise ValueError("Dish geometry yields zero/negative recess at inner radius (check R, a, fuel_ri).")
        if dz_inner >= 0.5 * pellet_height:
            raise ValueError(
                f"Dish too deep for pellet height: dz_inner={dz_inner:g} "
                f"must be < pellet_height/2={0.5*pellet_height:g}."
            )
        if chamfer_height <= 0.0 or chamfer_height >= 0.5 * pellet_height:
            raise ValueError("chamfer_height must satisfy 0 < chamfer_height < pellet_height/2.")

    else:
        raise ValueError(f"Unknown pellet_type: {pellet_type}")

    # -----------------------------
    # build mesh container
    # -----------------------------
    bm = Rod2DRZDiscreteBlockMesh(
        scale=scale,
        wedge_angle=wedge_angle,
        layout=layout,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        fuel_nr=fuel_nr,
        fuel_grading_nr=fuel_grading_nr,
        pellet_count=pellet_count,
        pellet_nz=pellet_nz,
        clad_nr=clad_nr,
        clad_nz=clad_nz,
        pellet_height=pellet_height,
        plenum_length=plenum_length,
        offset=offset,
        pellet_type=pellet_type,
        dish_radius_curvature=dish_radius_curvature,
        dish_outer_radius=dish_outer_radius,
        chamfer_height=chamfer_height,
        chamfer_start_radius=chamfer_start_radius,
        fuel_name=fuel_name,
        clad_name=clad_name,
    )

    # heights
    stack_height = pellet_count * pellet_height

    # -----------------------------
    # cladding: one continuous wedge
    # -----------------------------
    clad_block = None
    if has_clad:
        clad_height = stack_height + plenum_length
        clad_block = bm.create_wedge(
            clad_name,
            clad_ri,
            clad_ro,
            offset,
            offset + clad_height,
            wedgeAngle=wedge_angle,
            nr=clad_nr,
            nz=clad_nz,
        )

    # -----------------------------
    # fuel: pellets (each pellet has rings)
    # -----------------------------
    fuel_pellets: list[list] = []  # list[pellet][ring] -> Block

    if has_fuel:
        for i in range(pellet_count):
            z0 = offset + i * pellet_height
            z1 = z0 + pellet_height
            rings: list = []

            if pellet_type == "flat":
                b = bm.create_wedge(
                    fuel_name,
                    fuel_ri,
                    fuel_ro,
                    z0,
                    z1,
                    wedgeAngle=wedge_angle,
                    nr=fuel_nr,
                    nz=pellet_nz,
                    gradr=fuel_grading_nr,
                )
                rings.append(b)

            elif pellet_type == "dished":
                a = float(dish_outer_radius)
                R = float(dish_radius_curvature)

                nr0 = _nr_split(fuel_ri, a)
                nr1 = max(1, fuel_nr - nr0)

                b0 = bm.create_wedge_dished(
                    fuel_name,          # inner dished ring
                    fuel_ri,
                    a,
                    z0,
                    z1,
                    dishRadiusCurvature=R,
                    dishOuterRadius=a,
                    wedgeAngle=wedge_angle,
                    nr=nr0,
                    nz=pellet_nz,
                    gradr=fuel_grading_nr,
                )
                b1 = bm.create_wedge(
                    fuel_name,          # outer ring
                    a,
                    fuel_ro,
                    z0,
                    z1,
                    wedgeAngle=wedge_angle,
                    nr=nr1,
                    nz=pellet_nz,
                )
                rings.extend([b0, b1])

            elif pellet_type == "chamfered":
                rch = float(chamfer_start_radius)

                nr0 = _nr_split(fuel_ri, rch)
                nr1 = max(1, fuel_nr - nr0)

                b0 = bm.create_wedge(
                    fuel_name,
                    fuel_ri,
                    rch,
                    z0,
                    z1,
                    wedgeAngle=wedge_angle,
                    nr=nr0,
                    nz=pellet_nz,
                    gradr=fuel_grading_nr,
                )
                b1 = bm.create_wedge_chamfered(
                    fuel_name,
                    rch,
                    fuel_ro,
                    z0,
                    z1,
                    chamferHeight=float(chamfer_height),
                    wedgeAngle=wedge_angle,
                    nr=nr1,
                    nz=pellet_nz,
                )
                rings.extend([b0, b1])

            elif pellet_type == "dished_and_chamfered":
                a = float(dish_outer_radius)
                rch = float(chamfer_start_radius)
                R = float(dish_radius_curvature)

                nr0 = _nr_split(fuel_ri, a)
                nr1 = _nr_split(a, rch)
                nr2 = max(1, fuel_nr - (nr0 + nr1))

                b0 = bm.create_wedge_dished(
                    fuel_name,          # dish
                    fuel_ri,
                    a,
                    z0,
                    z1,
                    dishRadiusCurvature=R,
                    dishOuterRadius=a,
                    wedgeAngle=wedge_angle,
                    nr=nr0,
                    nz=pellet_nz,
                    gradr=fuel_grading_nr,
                )
                b1 = bm.create_wedge(
                    fuel_name,          # middle ring, sometimes called "land"
                    a,
                    rch,
                    z0,
                    z1,
                    wedgeAngle=wedge_angle,
                    nr=nr1,
                    nz=pellet_nz,
                )
                b2 = bm.create_wedge_chamfered(
                    fuel_name,          # chamfer ring
                    rch,
                    fuel_ro,
                    z0,
                    z1,
                    chamferHeight=float(chamfer_height),
                    wedgeAngle=wedge_angle,
                    nr=nr2,
                    nz=pellet_nz,
                )
                rings.extend([b0, b1, b2])

            fuel_pellets.append(rings)

            # --- internal patches for ring interfaces (per pellet) ---
            if len(rings) > 1:
                for k in range(len(rings) - 1):
                    f_out = mesh.Face(name=f"{fuel_name}Outer_{i}_{k}", boundaryType="patch")
                    f_out.add_sub_face(rings[k].rightFace())
                    bm.add_boundary(f_out)

                    f_in = mesh.Face(name=f"{fuel_name}Inner_{i}_{k+1}", boundaryType="patch")
                    f_in.add_sub_face(rings[k+1].leftFace())
                    bm.add_boundary(f_in)

            # --- internal patches for axial interfaces (pellet-to-pellet) ---
            if i > 0 and pellet_count > 1:
                f = mesh.Face(
                    name=f"{fuel_name}Bottom_{i}", 
                    boundaryType="regionCoupledOFFBEAT",
                    extraParameters={
                        "neighbourPatch": f"{fuel_name}Top_{i-1}",
                        "neighbourRegion": neighbour_region,
                        "owner": "false",
                        "updateAMI": "false",
                    },
                )
                for b in rings:
                    f.add_sub_face(b.bottomFace())
                bm.add_boundary(f)

            if i < (pellet_count - 1) and pellet_count > 1:
                f = mesh.Face(
                    name=f"{fuel_name}Top_{i}", 
                    boundaryType="regionCoupledOFFBEAT",
                    extraParameters={
                        "neighbourPatch": f"{fuel_name}Bottom_{i+1}",
                        "neighbourRegion": neighbour_region,
                        "owner": "true",
                        "updateAMI": "false",
                    },
                )
                for b in rings:
                    f.add_sub_face(b.topFace())
                bm.add_boundary(f)

    # -----------------------------
    # external wedge faces
    # -----------------------------
    if has_fuel:
        bc_fuel_front = mesh.Face(name="fuelFront", boundaryType="wedge")
        bc_fuel_back = mesh.Face(name="fuelBack", boundaryType="wedge")
        for rings in fuel_pellets:
            for b in rings:
                bc_fuel_front.add_sub_face(b.frontFace())
                bc_fuel_back.add_sub_face(b.backFace())
        bm.add_boundary(bc_fuel_front)
        bm.add_boundary(bc_fuel_back)

    if has_clad:
        bc_clad_front = mesh.Face(name="cladFront", boundaryType="wedge")
        bc_clad_front.add_sub_face(clad_block.frontFace())
        bm.add_boundary(bc_clad_front)

        bc_clad_back = mesh.Face(name="cladBack", boundaryType="wedge")
        bc_clad_back.add_sub_face(clad_block.backFace())
        bm.add_boundary(bc_clad_back)

    # -----------------------------
    # external radial boundaries
    # -----------------------------
    if has_fuel and fuel_ri > 0.0:
        bc_fuel_inner = mesh.Face(name="fuelInner")
        for rings in fuel_pellets:
            bc_fuel_inner.add_sub_face(rings[0].leftFace())
        bm.add_boundary(bc_fuel_inner)

    if has_fuel and has_clad:
        update_ami_str = "true" if updateAMI else "false"
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
        for rings in fuel_pellets:
            bc_fuel_outer.add_sub_face(rings[-1].rightFace())
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
        for rings in fuel_pellets:
            bc_fuel_outer.add_sub_face(rings[-1].rightFace())
        bm.add_boundary(bc_fuel_outer)

    elif has_clad and not has_fuel:
        bc_clad_inner = mesh.Face(name="cladInner")
        bc_clad_inner.add_sub_face(clad_block.leftFace())
        bm.add_boundary(bc_clad_inner)

        bc_clad_outer = mesh.Face(name="cladOuter")
        bc_clad_outer.add_sub_face(clad_block.rightFace())
        bm.add_boundary(bc_clad_outer)

    # -----------------------------
    # external caps (keep only global top/bottom)
    # -----------------------------
    if has_fuel:
        bc_fuel_bottom = mesh.Face(name="fuelBottom", boundaryType="patch")
        for b in fuel_pellets[0]:
            bc_fuel_bottom.add_sub_face(b.bottomFace())
        bm.add_boundary(bc_fuel_bottom)

        bc_fuel_top = mesh.Face(name="fuelTop", boundaryType="patch")
        for b in fuel_pellets[-1]:
            bc_fuel_top.add_sub_face(b.topFace())
        bm.add_boundary(bc_fuel_top)

    if has_clad:
        bc_clad_bottom = mesh.Face(name="cladBottom", boundaryType="patch")
        bc_clad_bottom.add_sub_face(clad_block.bottomFace())
        bm.add_boundary(bc_clad_bottom)

        bc_clad_top = mesh.Face(name="cladTop", boundaryType="patch")
        bc_clad_top.add_sub_face(clad_block.topFace())
        bm.add_boundary(bc_clad_top)

    # -----------------------------
    # deterministic mergePatchPairs (cheap)
    # -----------------------------
    if has_fuel and merge_radial_rings:
        for i, rings in enumerate(fuel_pellets):
            if len(rings) > 1:
                for k in range(len(rings) - 1):
                    bm.mergePatchPairsByName(
                        f"{fuel_name}Outer_{i}_{k}",
                        f"{fuel_name}Inner_{i}_{k+1}",
                    )
 
    return bm
