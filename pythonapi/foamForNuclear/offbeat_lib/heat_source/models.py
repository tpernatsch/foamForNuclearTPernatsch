# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

from ..profiles import axial_profile
from ..profiles import radial_profile

@offbeat_define
class HeatSource(OffbeatDict):
    """
    Base heat source model class. It neglects the heat source entirely, allowing
    for a simulation where the heat source field is not at all considered.

    This means that the `Q` field is not created.
    """
    TYPE: ClassVar[str] = 'none'

@offbeat_define
class Constant(HeatSource):
    """
    Heat source model class that allows the user to set the volumetric heat source
    field by providing or setting the corresponding field file `Q` in the starting
    time folder.

    If the `Q` file is present in the starting time folder, the internal field and boundary
    conditions are set according to its specifications. Otherwise, the heat source
    is set to 0 W/m3 in each cell, with zeroGradient BCs in all non-empty/-wedge patches.

    Note that when selecting this heat source model, the heat source field will remain
    the same as it was defined in the starting time folder (or as it was set by default
    by the code).
    """
    TYPE: ClassVar[str] = 'constant'

@offbeat_define
class ConstantLhgr(Constant):
    """
    Heat source model class that allows for a constant average **linear heat
    generation rate (lhgr)**.

    It is possible to apply a radial and an axial profile to the lhgr. Note
    that even if the lhgr is constant over time, the axial and radial profile might
    change over time, depending on the profile type chosen by the user.

    It is also possible to use the external reference dimensions option to apply the
    lhgr of the full-lenght rod a small segment or rodlet.


    Options
    -------
    lhgr : scalar
        The constant linear heat generation rate in W/m.
        (required: True)

    materials : list[word]
        A list of material (or cellZones names) where the heat source model applies.
        (required: True)

    axialProfile : axialProfile
        A sub-dictionary that specifies the type of axial profile.
        (required: True)

    radialProfile : radialProfile
        A sub-dictionary that specifies the type of radial profile.
        (required: True)

    useExternalReferenceDimensions : bool
        Activate reading reference dimensions for axial profile (`zMin` and `zMax`) from the
        dictionary.
        (default: False; required: False)

    zMin : scalar
        Reference axial coordinate of the bottom of the fuel region (only if
        `useReferenceDimensions` is true).
        (required: False)

    zMax : scalar
        Reference axial coordinate of the top of the fuel region (only if
        `useReferenceDimensions` is true).
        (required: False)

    volumeFraction : scalar
        Fraction of the reference rod present in the local mesh (only if
        `useReferenceDimensions` is true).
        (required: False)
    """
    TYPE: ClassVar[str] = 'constantLhgr'
    lhgr: float | int
    materials: list[str]
    axialProfile: axial_profile.AxialProfile = field(factory=axial_profile.Flat)
    radialProfile: radial_profile.RadialProfile = field(factory=radial_profile.Flat)
    useExternalReferenceDimensions: bool = False
    zMin: float | int | None = None
    zMax: float | int | None = None
    volumeFraction: float | int | None = None

@offbeat_define
class LaserHeatSource(Constant):
    """
    Heat source model class used to simulate the heating due to energy
    deposition from an impinging laser on one of the surfaces of the domain. The
    model initializes a `volScalarField` for the laser intensity, called I.


    Options
    -------
    alphaAbsorption : scalar
        The absorption coefficient of the material in $m^{-1}$.
        (required: True)

    beamDirection : vector
        A vector that specifies the direction of the laser beam.
        (required: True)

    useEffectiveAlpha : bool
        Whether to use the effective absorption coefficient.
        (required: True)

    materials : list[word]
        A list of material (or cellZones names) where the heat source model applies.
        (required: True)
    """
    TYPE: ClassVar[str] = 'laserHeatSource'
    alphaAbsorption: float | int
    beamDirection: list[int | float]
    useEffectiveAlpha: bool
    materials: list[str]

@offbeat_define
class TimeDependentVhgr(Constant):
    """
    Heat source model class that allows the user to set an average **volumetric
    heat generation rate (vhgr)**. The vhgr provided by the user can vary over time.

    This class is mainly developed for TRISO, and it currently doesn't support the
    addition of a radial and an axial profile of vhgr (as done in the lhgr heat
    source classes).


    Options
    -------
    timePoints : list[scalar]
        A list of time values at which the vhgr is provided. The time unit depends
        on the userTime selected by the user (seconds by default).
        (required: True)

    vhgr : list[scalar]
        A list of vhgr values in W/m3, one value per time-point.
        (required: True)

    timeInterpolationMethod : word
        Select the time interpolation method for time steps that fall in between the time
        points in `timePoints`.
        (default: 'linear'; required: True)

    materials : list[word]
        A list of material (or cellZones names) where the heat source model applies.
        (required: True)
    """
    TYPE: ClassVar[str] = 'timeDependentVhgr'
    timePoints: list[float | int]
    vhgr: list[float | int]
    timeInterpolationMethod: str = 'linear'
    materials: list[str]

@offbeat_define
class TimeDependentLhgr(ConstantLhgr):
    """
    Heat source model class that, similarly to the `constantLhgr` heatSource class,
    allows the user to set a average **linear heat generation rate (lhgr)**.

    The main difference is that the lhgr provided by the user can vary over time.

    It is possible to apply a radial and an axial profile to the lhgr. Note
    that even if the lhgr is constant over time, the axial and radial profile might
    change over time, depending on the profile type chosen by the user.

    It is also possible to use the external reference dimensions option to apply the
    lhgr of the full-lenght rod a small segment or rodlet.


    Options
    -------
    timePoints : list[scalar]
        A list of time values at which the lhgr is provided. The time unit depends
        on the userTime selected by the user (seconds by default).
        (required: True)

    lhgr : list[scalar]
        A list of lhgr values in W/m, one value per time-point.
        (required: True)

    timeInterpolationMethod : word
        Select the time interpolation method for time steps that fall in between the time
        points in `timePoints`.
        (default: 'linear'; required: True)

    materials : list[word]
        A list of material (or cellZones names) where the heat source model applies.
        (required: True)

    axialProfile : axialProfile
        A sub-dictionary that specifies the type of axial profile.
        (required: True)

    radialProfile : radialProfile
        A sub-dictionary that specifies the type of radial profile.
        (required: True)

    useExternalReferenceDimensions : bool
        Activate reading reference dimensions for axial profile (`zMin` and `zMax`) from the
        dictionary.
        (default: False; required: False)

    zMin : scalar
        Reference axial coordinate of the bottom of the fuel region (only if
        `useReferenceDimensions` is true).
        (required: False)

    zMax : scalar
        Reference axial coordinate of the top of the fuel region (only if
        `useReferenceDimensions` is true).
        (required: False)

    volumeFraction : scalar
        Fraction of the reference rod present in the local mesh (only if
        `useReferenceDimensions` is true).
        (required: False)
    """
    TYPE: ClassVar[str] = 'timeDependentLhgr'
    timePoints: list[float | int]
    lhgr: list[float | int]
    timeInterpolationMethod: str = 'linear'
    materials: list[str]
    axialProfile: axial_profile.AxialProfile = field(factory=axial_profile.Flat)
    radialProfile: radial_profile.RadialProfile = field(factory=radial_profile.Flat)
    useExternalReferenceDimensions: bool = False
    zMin: float | int | None = None
    zMax: float | int | None = None
    volumeFraction: float | int | None = None
