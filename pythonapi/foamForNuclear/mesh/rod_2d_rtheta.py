from __future__ import annotations

from typing import Literal

import attrs
import foamForNuclear.mesh as mesh


@attrs.define(init=False, slots=False, kw_only=True)
class Rod2DRThetaBlockMesh(mesh.BlockMesh):
    scale: float = 1.0

    fuel_ri: float
    fuel_ro: float
    fuel_nr: int
    fuel_nt: int

    clad_ri: float
    clad_ro: float
    clad_nr: int
    clad_nt: int

    z0: float = 0.0
    z1: float = 1.0
    nz: int = 1

    fuel_name: str = "fuel"
    clad_name: str = "cladding"

    x: float = 0.0
    y: float = 0.0

    def __init__(
        self,
        *,
        scale: float = 1.0,
        fuel_ri: float,
        fuel_ro: float,
        fuel_nr: int,
        fuel_nt: int,
        clad_ri: float,
        clad_ro: float,
        clad_nr: int,
        clad_nt: int,
        z0: float = 0.0,
        z1: float = 1.0,
        nz: int = 1,
        fuel_name: str = "fuel",
        clad_name: str = "cladding",
        x: float = 0.0,
        y: float = 0.0,
    ):
        super().__init__(scale=scale)

        self.scale = scale
        self.fuel_ri = fuel_ri
        self.fuel_ro = fuel_ro
        self.fuel_nr = fuel_nr
        self.fuel_nt = fuel_nt

        self.clad_ri = clad_ri
        self.clad_ro = clad_ro
        self.clad_nr = clad_nr
        self.clad_nt = clad_nt

        self.z0 = z0
        self.z1 = z1
        self.nz = nz

        self.fuel_name = fuel_name
        self.clad_name = clad_name

        self.x = x
        self.y = y


def rod_2d_rtheta(
    *,
    scale: float = 1.0,
    # fuel ring
    fuel_ri: float,
    fuel_ro: float,
    fuel_nr: int,
    fuel_nt: int,
    # clad ring
    clad_ri: float,
    clad_ro: float,
    clad_nr: int,
    clad_nt: int,
    fuel_name: str = "fuel",
    clad_name: str = "cladding",
    # 2D extrusion in z
    z0: float = 0.0,
    z1: float = 1.0,
    nz: int = 1,
    # clad center
    x: float = 0.0,
    y: float = 0.0,
    # NEW: eccentricity of fuel center wrt clad center
    eccentricity: float = 0.0,
    ecc_dir: Literal["x", "y"] = "x",
    neighbour_region: str = "region0",
    updateAMI: bool = False,
) -> Rod2DRThetaBlockMesh:
    """
    2D r-theta (x-y plane) concentric fuel+clad rings.
    top/bottom are 'empty' and named exactly 'top'/'bottom'.
    fuelOuter <-> cladInner are regionCoupledOFFBEAT patches.

    Eccentricity:
      - shifts the fuel ring center by `eccentricity` along `ecc_dir`
        relative to the clad center (x, y).
      - must satisfy: eccentricity <= clad_ri - fuel_ro  (no overlap)
    """

    # --- checks ---
    if nz != 1:
        raise ValueError("For a 2D r-theta mesh, nz must be 1 (top/bottom are empty).")
    if z1 == z0:
        raise ValueError("z1 must be different from z0 (non-zero extrusion thickness).")

    if not (fuel_ri > 0.0 and fuel_ro > fuel_ri):
        raise ValueError("Fuel radii must satisfy: fuel_ri > 0 and fuel_ro > fuel_ri.")
    if not (clad_ri > 0.0 and clad_ro > clad_ri):
        raise ValueError("Cladding radii must satisfy: clad_ri > 0 and clad_ro > clad_ri.")
    if clad_ri < fuel_ro:
        raise ValueError("Expected clad_ri >= fuel_ro (gap >= 0, no overlap).")

    if fuel_nr <= 0 or fuel_nt <= 0 or clad_nr <= 0 or clad_nt <= 0:
        raise ValueError("nr/nt must be positive integers.")

    if eccentricity < 0.0:
        raise ValueError("eccentricity must be >= 0.")
    gap = clad_ri - fuel_ro
    if eccentricity > gap:
        raise ValueError(
            f"eccentricity={eccentricity:g} exceeds available gap={gap:g} "
            "(need eccentricity <= clad_ri - fuel_ro)."
        )

    bm = Rod2DRThetaBlockMesh(
        scale=scale,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        fuel_nr=fuel_nr,
        fuel_nt=fuel_nt,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        clad_nr=clad_nr,
        clad_nt=clad_nt,
        z0=z0,
        z1=z1,
        nz=nz,
        fuel_name=fuel_name,
        clad_name=clad_name,
        x=x,
        y=y,
    )

    # --- centers ---
    x_clad, y_clad = x, y
    if ecc_dir == "x":
        x_fuel, y_fuel = x + eccentricity, y
    else:
        x_fuel, y_fuel = x, y + eccentricity

    # --- create two rings (each = 4 blocks) ---
    fuel_blocks = bm.create_ring_along_z(
        name=fuel_name,
        innerRadius=fuel_ri,
        outerRadius=fuel_ro,
        lowZ=z0,
        highZ=z1,
        nr=fuel_nr,
        nt=fuel_nt,
        nz=nz,
        x=x_fuel,
        y=y_fuel,
        isAddBoundaryConditions=False,
    )
    clad_blocks = bm.create_ring_along_z(
        name=clad_name,
        innerRadius=clad_ri,
        outerRadius=clad_ro,
        lowZ=z0,
        highZ=z1,
        nr=clad_nr,
        nt=clad_nt,
        nz=nz,
        x=x_clad,
        y=y_clad,
        isAddBoundaryConditions=False,
    )

    fuel_front, fuel_right, fuel_back, fuel_left = fuel_blocks
    clad_front, clad_right, clad_back, clad_left = clad_blocks

    # --- top/bottom as EMPTY ---
    top = mesh.Face(name="top", boundaryType="empty")
    bottom = mesh.Face(name="bottom", boundaryType="empty")
    for b in [*fuel_blocks, *clad_blocks]:
        top.add_sub_face(b.topFace())
        bottom.add_sub_face(b.bottomFace())
    bm.add_boundary(top)
    bm.add_boundary(bottom)

    # --- coupled interface fuelOuter <-> cladInner ---
    update_ami_str = "true" if updateAMI else "false"

    fuelOuter = mesh.Face(
        name="fuelOuter",
        boundaryType="regionCoupledOFFBEAT",
        extraParameters={
            "neighbourPatch": "cladInner",
            "neighbourRegion": neighbour_region,
            "owner": "true",
            "updateAMI": f"{update_ami_str}",
        },
    )
    fuelOuter.add_sub_face(fuel_front.frontFace())
    fuelOuter.add_sub_face(fuel_left.leftFace())
    fuelOuter.add_sub_face(fuel_right.rightFace())
    fuelOuter.add_sub_face(fuel_back.backFace())
    bm.add_boundary(fuelOuter)

    cladInner = mesh.Face(
        name="cladInner",
        boundaryType="regionCoupledOFFBEAT",
        extraParameters={
            "neighbourPatch": "fuelOuter",
            "neighbourRegion": neighbour_region,
            "owner": "false",
            "updateAMI": f"{update_ami_str}",
        },
    )
    cladInner.add_sub_face(clad_front.backFace())
    cladInner.add_sub_face(clad_left.rightFace())
    cladInner.add_sub_face(clad_right.leftFace())
    cladInner.add_sub_face(clad_back.frontFace())
    bm.add_boundary(cladInner)

    # --- outer boundary of cladding ---
    cladOuter = mesh.Face(name="cladOuter")
    cladOuter.add_sub_face(clad_front.frontFace())
    cladOuter.add_sub_face(clad_left.leftFace())
    cladOuter.add_sub_face(clad_right.rightFace())
    cladOuter.add_sub_face(clad_back.backFace())
    bm.add_boundary(cladOuter)

    # --- inner boundary of fuel (central hole) ---
    fuelInner = mesh.Face(name="fuelInner")
    fuelInner.add_sub_face(fuel_front.backFace())
    fuelInner.add_sub_face(fuel_left.rightFace())
    fuelInner.add_sub_face(fuel_right.leftFace())
    fuelInner.add_sub_face(fuel_back.frontFace())
    bm.add_boundary(fuelInner)

    return bm
