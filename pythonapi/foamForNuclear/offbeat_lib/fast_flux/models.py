# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

from ..profiles import axial_profile

@offbeat_define
class FastFlux(OffbeatDict):
    """
    Base fast flux model class. It disables entirely the fast flux and fast
    fluence handling in OFFBEAT.

    When selected, the fast flux and fast fluence fields (`fastFlux`, `fastFluence`)
    are **not created**, meaning these quantities are not available to other models.

    This is the default choice when no fast flux model is explicitly selected.
    """
    TYPE: ClassVar[str] = 'none'

@offbeat_define
class Constant(FastFlux):
    """
    Fast flux model that allows prescribing a fixed fast neutron flux
    field in OFFBEAT via the `fastFlux` field file in the starting folder.

    When this model is selected, the fast flux (`fastFlux`) and fast fluence
    (`fastFluence`) fields are read from the starting time folder if present.
    If the corresponding files are not found, both fields are initialized to zero:
    - `fastFlux` is set to 0 n/cm²/s
    - `fastFluence` is set to 0 n/cm²

    In this case, zeroGradient boundary conditions are applied to all non-empty
    and non-wedge patches.

    The fast flux field remains constant throughout the simulation, i.e. it is
    not modified during time advancement.

    Even though the fast flux is constant, the fast fluence is updated at every
    time step by integrating the fast flux over time.
    """
    TYPE: ClassVar[str] = 'constant'

@offbeat_define
class TimeDependentAxialProfile(Constant):
    """
    Fast flux model that allows prescribing a **time-dependent rod-average
    fast neutron flux**, optionally combined with a **time- and space-dependent
    axial profile**.

    The fast fluence evolves in time according to the local value of the fast flux.
    Even if the rod-average fast flux is constant, the axial profile may vary in
    time depending on the selected profile model.


    Options
    -------
    timePoints : list[scalar]
        List of time points at which the rod-average fast flux values are provided.
        The time unit depends on the `userTime` setting (seconds by default).
        (required: True)

    fastFlux : list[scalar]
        List of rod-average fast flux values (n/cm²/s), one per entry in
        `timePoints`.
        (required: True)

    timeInterpolationMethod : word
        Time interpolation method used between entries in `timePoints`.
        Typical options are `linear` or `step`.
        (default: 'linear'; required: False)

    materials : list[word]
        List of materials or cellZones where the fast flux model is applied.
        If omitted, the model is applied everywhere.
        (required: False)

    axialProfile : axialProfile
        Sub-dictionary defining the axial profile model.

        The profile may be time-independent (e.g. `flat`) or time-dependent
        (e.g. `timeDependentTabulated`). Refer to the `axialProfile` class documentation
        for available models and required parameters.
        (required: True)
    """
    TYPE: ClassVar[str] = 'timeDependentAxialProfile'
    timePoints: list[float | int]
    fastFlux: list[float | int]
    timeInterpolationMethod: str = 'linear'
    materials: list[str] | None = None
    axialProfile: axial_profile.AxialProfile = field(factory=axial_profile.Flat)
