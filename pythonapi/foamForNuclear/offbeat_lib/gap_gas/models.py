# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

from foamForNuclear.timeProfile import OffbeatTimeProfile

@offbeat_define
class GapGasModel(OffbeatDict):
    """
    Base gap gas model class. It neglects the presence of any gap gas.
    No gas pressure, composition, or thermal transport through the gap is modeled.
    """
    TYPE: ClassVar[str] = 'none'

@offbeat_define
class Frapcon(GapGasModel):
    """
    Gap gas model that computes gap gas composition, pressure, temperature,
    and volume using correlations derived from the FRAPCON 4.0 manual.

    The gap gas state at the initial time (presure and composition) must
    be separately provided via the `gapGas` file.

    The model assumes:
    - ideal gas behavior,
    - multiple free-volume regions (gap, plena, central hole, cracks),
    - instantaneous pressure equalization over the entire free volume.


    Options
    -------
    gapPatches : list[word]
        Patches defining the fuel–cladding gap (e.g. `fuelOuter cladInner`).
        (required: True)

    holePatches : list[word]
        Patches defining internal fuel holes.
        (default: []; required: False)

    topFuelPatches : list[word]
        Patches at the top of the fuel stack.
        (required: False)

    bottomFuelPatches : list[word]
        Patches at the bottom of the fuel stack.
        (required: False)

    gapVolumeOffset : scalar
        Offset applied to the initial gap volume.
        (default: 0.0; required: False)

    gasReserveVolume : scalar
        Additional reserve gas volume used in pressure calculation.
        (default: 0.0; required: False)

    gasReserveTemperature : scalar
        Temperature of the reserve gas (K).
        (default: 290; required: False)

    gasReserveTemperatureList : table
        Time-dependent table for reserve gas temperature.
        If present, it overrides `gasReserveTemperature`.
        (required: False)

    includeCentralHole : bool
        Include fuel central hole volume in the free volume calculation.
        (default: True; required: False)

    includeDishes : bool
        Include dish geometry in the free volume calculation.
        (default: True; required: False)
    """
    TYPE: ClassVar[str] = 'Frapcon'
    gapPatches: list[str]
    holePatches: list[str] = field(factory=list)
    topFuelPatches: list[str] | None = None
    bottomFuelPatches: list[str] | None = None
    gapVolumeOffset: float | int = 0.0
    gasReserveVolume: float | int = 0.0
    gasReserveTemperature: float | int = 290.0
    gasReserveTemperatureList: OffbeatTimeProfile | None = None
    includeCentralHole: bool = True
    includeDishes: bool = True

@offbeat_define
class TimeTabulated(GapGasModel):
    """
    Gap gas model class that prescribes **time-dependent gap gas
    composition and pressure** from externally provided data.

    This model is intended for coupling with external fuel performance codes
    or for replaying precomputed gap gas histories.

    A snapshot of the gap gas state at the initial time must still be provided
    via the `gapGas` file.


    Options
    -------
    timePoints : list[scalar]
        Time points corresponding to the tabulated gap gas data.
        (required: True)

    gapPressures : list[scalar]
        Gap gas pressures corresponding to `timePoints` (Pa).
        (required: True)

    Ar : list[scalar]
        Argon mass fraction history.
        (required: True)

    He : list[scalar]
        Helium mass fraction history.
        (required: True)

    Kr : list[scalar]
        Krypton mass fraction history.
        (required: True)

    Ne : list[scalar]
        Neon mass fraction history.
        (required: True)

    Rn : list[scalar]
        Radon mass fraction history.
        (required: True)

    Xe : list[scalar]
        Xenon mass fraction history.
        (required: True)

    timeInterpolationMethod : word
        Time interpolation scheme between tabulated values.
        (default: 'linear'; required: False)
    """
    TYPE: ClassVar[str] = 'timeTabulated'
    timePoints: list[float | int]
    gapPressures: list[float | int]
    Ar: list[float | int]
    He: list[float | int]
    Kr: list[float | int]
    Ne: list[float | int]
    Rn: list[float | int]
    Xe: list[float | int]
    timeInterpolationMethod: str = 'linear'

@offbeat_define
class TRISO(GapGasModel):
    """
    Gap gas model class that computes gas composition, volume, temperature, and
    pressure for TRISO fuel configurations.

    The gap gas state at the initial time (presure and composition) must
    be separately provided via the `gapGas` file.

    The model is similar to the `Frapcon` gap gas model in terms of assumptions
    and gap volume evaluation strategy (ideal gas, free volume decomposition,
    and volume estimation from bounding surfaces). See `Frapcon` documentation
    for details on the gap volume algorithm.


    Options
    -------
    heatSourceName : word
        Name of the volumetric heat source field used to infer fission rate / CO
        production-related quantities.
        (default: 'Q'; required: False)

    COProductionModel : word
        CO production model selection: `Proksch` (default) or `GA`.
        (default: 'Proksch'; required: False)

    model3D : bool
        If true, interpret the modeled configuration as 3D.
        (default: False; required: False)

    bufferPorosity : scalar
        Porosity of the buffer layer.
        (required: True)

    bufferZones : list[word]
        Names of the cellZones corresponding to the buffer material.
        (required: True)

    gapPatches : list[word]
        Patches defining the TRISO gap (e.g. `bufferOuter ipycInner`).
        (required: True)

    gapVolumeOffset : scalar
        Offset applied to the initial gap volume.
        (default: 0.0; required: False)

    gasReserveVolume : scalar
        Initial volume of reserve gas for pressure calculation.
        (default: 0.0; required: False)

    gasReserveTemperature : scalar
        Initial temperature of reserve gas (K).
        (default: 290; required: False)

    gasReserveTemperatureList : table
        Optional time-dependent table for reserve gas temperature.
        If present, it overrides `gasReserveTemperature`.

        Expected structure:
        - `type`: table
        - `values`: list of (time, temperature)
        - optional `outOfBounds` (default: clamp)
        - optional `interpolationScheme` (e.g. linear)
        (required: False)

    includeCentralHole : bool
        Toggle inclusion of central hole volume in the free volume calculation.
        (default: True; required: False)

    includeDishes : bool
        Toggle inclusion of dish geometry in the free volume calculation.
        (default: True; required: False)
    """
    TYPE: ClassVar[str] = 'TRISO'
    heatSourceName: str = 'Q'
    COProductionModel: str = 'Proksch'
    model3D: bool = False
    bufferPorosity: float | int
    bufferZones: list[str]
    gapPatches: list[str]
    gapVolumeOffset: float | int = 0.0
    gasReserveVolume: float | int = 0.0
    gasReserveTemperature: float | int = 290.0
    gasReserveTemperatureList: OffbeatTimeProfile | None = None
    includeCentralHole: bool = True
    includeDishes: bool = True
