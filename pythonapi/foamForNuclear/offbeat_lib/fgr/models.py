# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

from foamForNuclear.common import SciantixDict

@offbeat_define
class FissionGasRelease(OffbeatDict):
    """
    Base fission gas release (FGR) model class. It disables fission gas production and
    release, and therefore also disables any FGR-related gaseous swelling.

    With this selection, **no additional fission-gas fields** (e.g. concentrations,
    release fractions) are created or evolved.
    """
    TYPE: ClassVar[str] = 'none'

@offbeat_define
class SCIANTIX(FissionGasRelease):
    """
    Fission gas release (FGR) model that couples OFFBEAT with the 0-D
    SCIANTIX code developed at Politecnico di Milano.


    Options
    -------
    relax : scalar
        Relaxation factor applied to fission gas release and gaseous swelling
        quantities when coupling SCIANTIX with OFFBEAT.
        (default: 1; required: False)

    nFrequency : int
        Number of outer iterations to skip between successive SCIANTIX executions.
        Increasing this value reduces computational cost but may affect convergence.
        (default: 1; required: False)

    addToRegistry : bool
        If enabled, SCIANTIX variables are registered as volumetric OpenFOAM fields
        (e.g. `volScalarField`) and written to disk, allowing visualization and
        probing in ParaView.
        (default: False; required: False)

    SCIANTIX : dictionary
        Sub-dictionary defining SCIANTIX input parameters.

        If present, this dictionary overrides the external `input_settings.txt`
        file. Parameters correspond directly to SCIANTIX model switches and
        numerical options (e.g. grain growth, diffusion models, bubble evolution,
        solver selection).
        (required: False)
    """
    TYPE: ClassVar[str] = 'SCIANTIX'
    relax: float | int = 1.0
    nFrequency: int = 1
    addToRegistry: bool = False
    SCIANTIX: SciantixDict = field(factory=SciantixDict)

@offbeat_define
class SCIANTIXRIA(FissionGasRelease):
    """
    Fission gas release (FGR) model that is intended for **RIA or
    transient simulations following a base irradiation phase**.

    The expected workflow is:
    1. Run the base irradiation using the `SCIANTIX` FGR model.
    2. Restart the simulation from a written time step.
    3. Switch the FGR model to `SCIANTIXRIA`.


    Options
    -------
    relax : scalar
        Relaxation factor applied to released gas fractions and associated
        swelling terms. Useful for improving numerical stability in strongly
        coupled transient simulations.
        (default: 1; required: False)

    releaseHBS : bool
        Activate fission gas release associated with High Burnup Structure (HBS)
        formation.

        When enabled, both burnup- and temperature-based thresholds are evaluated.
        (default: True; required: False)

    buReleaseThresholdHBS : scalar
        Burnup threshold for HBS-related gas release, expressed in MWd/t.
        (default: 80000; required: False)

    temperatureReleaseThresholdHBS : scalar
        Temperature threshold for HBS-related gas release, expressed in Kelvin.
        (default: 1000; required: False)

    damageReleaseThreshold : scalar
        Damage threshold triggering inter-granular gas release.
        The damage field is typically computed by the selected `damageModel`
        in the constitutive law.
        (default: 0.85; required: False)
    """
    TYPE: ClassVar[str] = 'SCIANTIXRIA'
    relax: float | int = 1.0
    releaseHBS: bool = True
    buReleaseThresholdHBS: float | int = 80000.0
    temperatureReleaseThresholdHBS: float | int = 1000.0
    damageReleaseThreshold: float | int = 0.85
