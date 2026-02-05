# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import offbeat_define
from foamForNuclear.common import OffbeatDict

@offbeat_define
class MultiMaterialInterface(OffbeatDict):
    """
    Base class for multi-material correction models.
    """
    TYPE: ClassVar[str] = 'multiMaterialInterface'

@offbeat_define
class Uniform(MultiMaterialInterface):
    """
    Multi-material correction class that applies a uniform multi-material correction.


    Options
    -------
    defaultWeights : scalar
        Fraction of the full correction applied (between 0 and 1).
        Applying full correction everywhere may negatively impact convergence.
        A value between 0.9 and 1 is recommended (1 if possible).
        (required: True)

    defaultWeightsGrad : scalar
        Fraction of the correction applied to the gradient contribution of
        `sigmaExp`. By default equal to `defaultWeights`.
        (required: False)
    """
    TYPE: ClassVar[str] = 'uniform'
    defaultWeights: float | int
    defaultWeightsGrad: float | int | None = None

@offbeat_define
class CellZone(Uniform):
    """
    Multi-material correction class where two sets of weights can be used:

    - `defaultWeights` for internal faces within a cellZone
    - `interfaceWeights` for interface faces between different cellZones


    Options
    -------
    defaultWeights : scalar
        Fraction of the correction applied on internal faces.
        Recommended between 0.9 and 1.
        (required: True)

    defaultWeightsGrad : scalar
        Fraction of the gradient correction.
        By default equal to `defaultWeights`.
        (required: False)

    interfaceWeights : scalar
        Fraction of the correction applied at interfaces between cellZones.
        Typically safe to set to 1.0.
        (required: True)
    """
    TYPE: ClassVar[str] = 'cellZone'
    defaultWeights: float | int
    defaultWeightsGrad: float | int | None = None
    interfaceWeights: float | int

@offbeat_define
class FaceSet(Uniform):
    """
    Multi-material correction class that uses a user-defined `faceSet`
    to identify interface faces.


    Options
    -------
    defaultWeights : scalar
        Fraction of the correction applied on internal faces.
        (required: True)

    defaultWeightsGrad : scalar
        Fraction of the gradient correction.
        (required: False)

    interfaceWeights : scalar
        Fraction of the correction applied at interface faces.
        (required: True)

    faceSetName : word
        Name of the faceSet identifying interface faces.
        (required: True)
    """
    TYPE: ClassVar[str] = 'faceSet'
    defaultWeights: float | int
    defaultWeightsGrad: float | int | None = None
    interfaceWeights: float | int
    faceSetName: str

@offbeat_define
class TopoSetSource(Uniform):
    """
    Multi-material correction class where the interface faces are defined from a
    `topoSetSource`. The correction is applied as follows:

    - first, the uniform correction is applied everywhere (via `defaultWeights`
    and optionally `defaultWeightsGrad`);
    - then, faces selected by the `source` entry are overwritten and forced to use
    `interfaceWeights`.


    Options
    -------
    defaultWeights : scalar
        Base (uniform) correction weight applied to internal faces not belonging
        to the selected interface.
        (required: True)

    defaultWeightsGrad : scalar
        Base correction weight for the gradient contribution. Defaults to
        `defaultWeights`.
        (required: False)

    interfaceWeights : scalar
        Correction weight applied to faces selected by the `source` topoSetSource.
        It overwrites both weightsEps and weightsGrad for those faces.
        (required: True)

    source : topoSetSource
        Definition of the interface faces using OpenFOAM `topoSetSource` syntax.
        The first token is the topoSetSource type, followed by its parameters.
        (required: True)
    """
    TYPE: ClassVar[str] = 'topoSetSource'
    defaultWeights: float | int
    defaultWeightsGrad: float | int | None = None
    interfaceWeights: float | int
    source: TopoSetSource

@offbeat_define
class Uniform2D(Uniform):
    """
    Multi-material correction class that applies a uniform multi-material correction
    **only in the x–y plane**.

    It was designed for 2D r–z simulations with a coarse axial discretization.
    In such cases, the full 3D correction may cause convergence issues due to
    insufficient axial resolution.

    As the axial mesh is refined, the standard `uniform` correction can be used
    safely.


    Options
    -------
    defaultWeights : scalar
        Fraction of the full correction applied (between 0 and 1).
        Applying full correction everywhere may negatively impact convergence.
        A value between 0.9 and 1 is recommended (1 if possible).
        (required: True)

    defaultWeightsGrad : scalar
        Fraction of the correction applied to the gradient contribution of
        `sigmaExp`. By default equal to `defaultWeights`.
        (required: False)
    """
    TYPE: ClassVar[str] = 'uniform2D'
    defaultWeights: float | int
    defaultWeightsGrad: float | int | None = None

@offbeat_define
class UniformDirectional(Uniform):
    """
    Multi-material correction class where the interfaces between materials are defined
    using directionally-dependent correction weights.

    The correction may optionally be defined in a rotated coordinate system,
    allowing the user to specify principal material directions.


    Options
    -------
    rotated : bool
        If true, face normals are transformed to a rotated/local coordinate system
        before applying directional weights.
        (required: True)

    rotation : dict
        Coordinate rotation definition (required if `rotated` is true). Supports
        standard OpenFOAM `coordinateRotation` types; may be spatially varying.
        (required: False)

    defaultWeights : vector
        Directional weights applied to the strain-related correction term.
        The per-face weight is computed as `| nf ⊙ defaultWeights |`, where `nf` is
        the unit face normal (in rotated coordinates if `rotated` is true).
        (required: True)

    defaultWeightsGrad : vector
        Directional weights applied to the gradient contribution in `sigmaExp`.
        If not provided, it defaults to `defaultWeights`.
        (required: False)
    """
    TYPE: ClassVar[str] = 'uniformDirectional'
    rotated: bool
    rotation: dict[str, Any] | None = None
    defaultWeights: list[int | float]
    defaultWeightsGrad: list[int | float] | None = None
