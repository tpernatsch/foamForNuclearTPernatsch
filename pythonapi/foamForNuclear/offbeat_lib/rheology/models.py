from __future__ import annotations

from foamForNuclear._attrs_tools import ffn_define
from typing import ClassVar
from foamForNuclear.common import FoamForNuclearDict, Table
from attrs import field

from foamForNuclear.timeProfile import OffbeatTimeProfile


@ffn_define
class Rheology(FoamForNuclearDict):
    TYPE: ClassVar[str] = "NOT IMPLEMENTED"


@ffn_define
class ByMaterial(Rheology):
    TYPE: ClassVar[str] = "byMaterial"
    thermalExpansion: bool = True


# Alias class
class Standard(ByMaterial):
    """User-facing name for API."""
    pass

ByMaterial = Standard

@ffn_define
class PlaneStress(Standard):
    planeStress: bool = True


def coolant_pressure_profile(value):
    if value is None:
        return None

    if isinstance(value, OffbeatTimeProfile):
        return value

    # Allow (data, outOfBounds)
    if isinstance(value, tuple) and len(value) == 2:
        data, outOfBounds, interpolation = value
        return OffbeatTimeProfile(
            type="table",
            values=data,
            outOfBounds=outOfBounds,
            interpolationScheme=interpolation
        )

    # Allow dict-style spec
    if isinstance(value, dict):
        return OffbeatTimeProfile(
            type="table",
            **value,
        )

    # Default: list or Table
    return OffbeatTimeProfile(
        type="table",
        values=value,
    )


@ffn_define
class ModifiedPlaneStrain(Standard):
    modifiedPlaneStrain: bool = True
    precisionSpring: float | int = 1
    springModulus: float | int = 3500.0
    coolantPressureList: list | Table | OffbeatTimeProfile | None = field(
        default=None,
        converter=coolant_pressure_profile,
    )    
    penaltyFriction: float | int = 1.0
    frictionCoefficient: float | int = 0.0
    frictionSaturationPressure: float | int = 1e7
    relax: float | int = 0.001