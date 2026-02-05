# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class FailureModel(OffbeatDict):
    """
    Mother class for materials failure criteria.

    The base model has type `none`. When selected, failure is not modeled.


    Options
    -------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)
    """
    TYPE: ClassVar[str] = 'none'
    stopIfFailed: bool = False

Failure = FailureModel  # alias

@offbeat_define
class Combined(Failure):
    """
    Criterion used to combine multiple failure criteria. Failure is supposed to
    occur when one of the criteria declare failure.

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    criteria : list[word]
        List of failure criteria to combine. Failure is triggered if any of the
        listed criteria declares failure.
        (required: True)
    """
    TYPE: ClassVar[str] = 'combined'
    criteria: list[str]

@offbeat_define
class UO2Matpro(Failure):
    """
    Model to check failure of UO2 based on a limit melting temperature from
    Matpro v11.

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    Tmelt : scalar
        Melting temperature threshold [K]. If the local fuel temperature exceeds
        this value, the material is considered failed.
        (default: 3113.0; required: False)
    """
    TYPE: ClassVar[str] = 'UO2Matpro'
    Tmelt: float | int = 3113.0

@offbeat_define
class UPuO2Magni2020(Failure):
    """
    Model to check failure of UPuO2 based on a limit melting temperature from
    "Modelling and assessment of thermal conductivity and melting
    behaviour of MOX fuel for fast reactor applications - A.Magni et al.".
    (https://doi.org/10.1016/j.jnucmat.2020.152410).

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    burnupName : word
        Name of the burnup field to read from the mesh registry.
        (default: 'Bu'; required: False)

    Tmelt : scalar
        Reference UO2 melting temperature parameter [K] used by the model.
        (default: 3147.0; required: False)
    """
    TYPE: ClassVar[str] = 'UPuO2Magni2020'
    burnupName: str = 'Bu'
    Tmelt: float | int = 3147.0

@offbeat_define
class ZircaloySed(Failure):
    """
    Clad failure criterion during the PCMI phase of a RIA transient based on the strain
    energy density (SED).
    The complete description of this model is available at :
    https://www.epri.com/research/products/1021036
    If the average SED on an axial slice goes over a critical value, the clad is considered
    to have failed.
    The SED field corresponds to the strain energy density.
    The CFP field corresponds to the SED field over the critical value. Thus if CFP > 1,
    the clad is considered failed.

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    CH : scalar
        Hydrogen content in wppm. The critical value is dependent on the hydrogen content.
        (required: True)

    crackLength : scalar
        Crack length used by the model [m].
        (required: True)
    """
    TYPE: ClassVar[str] = 'ZircaloySed'
    CH: float | int
    crackLength: float | int

@offbeat_define
class ZircaloyOverstrainBison(Failure):
    """
    Failure criterion based on an overstrain limit for Zircaloy cladding, inspired
    by the implementation used in BISON.

    Failure is triggered when the equivalent plastic strain exceeds a critical
    threshold.

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    patchNames : list[word]
        List of patch names on which the criterion is evaluated.
        (required: True)

    hoopStrainLimit : scalar
        Critical equivalent plastic strain above which the cladding is considered
        failed.
        (default: 0.4; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyOverstrainBison'
    patchNames: list[str]
    hoopStrainLimit: float | int = 0.4

@offbeat_define
class ZircaloyOverstressBison(Failure):
    """
    Burst overstress criterion for Zircaloy-4 derived from Bison code. NOTE:
    oxidation contribution still needs to be included.

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    patchNames : list[word]
        List of patch names where the criterion applies.
        (required: True)
    """
    TYPE: ClassVar[str] = 'ZircaloyOverstressBison'
    patchNames: list[str]

@offbeat_define
class ZircaloyPlasticInstabilityBison(Failure):
    """
    Plastic instability failure criterion for Zircaloy cladding inspired by BISON.

    The model evaluates plastic instability on the specified patches and uses a
    burst strain-rate threshold as an input parameter.

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    patchNames : list[word]
        List of patch names on which the criterion is evaluated.
        (required: True)

    burstStrainRate : scalar
        Burst strain-rate threshold used by the model.
        (default: 0.0278; required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyPlasticInstabilityBison'
    patchNames: list[str]
    burstStrainRate: float | int = 0.0278

@offbeat_define
class ZircaloyRiaJernkvistModified(Failure):
    """
    Modified Jernkvist failure criterion for Zircaloy cladding during RIA
    transients.

    Options (inherited)
    -------------------
    stopIfFailed : bool
        If true, the simulation is stopped when the failure criterion is met.
        (default: False; required: False)

    Options
    -------
    phi : scalar
        Fast neutron fluence in 1e25 n/m2
        (required: True)

    meanHydrogenConcentration : scalar
        Mean hydrogen concentration in the cladding (wppm).
        (required: True)

    rimThicknessReductionFactor : bool
        If true, activates a rim-thickness reduction factor and requires additional
        rim-related inputs (`rimThickness`, `cladOuterRadius`, `rimHydrogenConcentration`).
        (default: False; required: False)

    rimThickness : scalar
        Rim thickness (in um) used when `rimThicknessReductionFactor` is enabled.
        (required: False)

    cladOuterRadius : scalar
        Cladding outer radius (in m) used when `rimThicknessReductionFactor` is enabled.
        (required: False)

    rimHydrogenConcentration : scalar
        Hydrogen concentration in the rim region (wppm), used when
        `rimThicknessReductionFactor` is enabled.
        (required: False)
    """
    TYPE: ClassVar[str] = 'ZircaloyRiaJernkvistModified'
    phi: float | int
    meanHydrogenConcentration: float | int
    rimThicknessReductionFactor: bool = False
    rimThickness: float | int | None = None
    cladOuterRadius: float | int | None = None
    rimHydrogenConcentration: float | int | None = None
