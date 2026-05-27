# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict

from . import multi_material

@ffn_define
class MechanicsSubSolver(FoamForNuclearDict):
    """
    Base mechanics solver class. It keeps the displacement field `D` unchanged
    from its initial configuration throughout the simulation.

    While the `constant` solver keeps the `D` field fixed by default, it is possible
    to alter the `D` field during the simulation by leveraging advanced OpenFOAM
    features like `#codeStream#`, `funkySetFields`, or specific function objects.


    Options
    -------
    forceSummary : bool
        Print a summary of boundary-forces at each iteration.
        (default: False; required: False)

    cylindricalStress : bool
        Creates cylindrical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    sphericalStress : bool
        Creates spherical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    RhieChowCorrection : bool
        Apply a Rhie–Chow–like correction to stabilize the momentum equation.
        (default: True; required: False)

    RhieChowScaleFactor : scalar
        Scaling factor for the Rhie–Chow correction.
        (default: 1; required: False)

    multiMaterialCorrection : dict
        Optional sub-dictionary activating a custom discretization for
        multi-material interfaces.
        (required: False)
    """
    TYPE: ClassVar[str] = 'constant'
    forceSummary: bool = False
    cylindricalStress: bool = False
    sphericalStress: bool = False
    RhieChowCorrection: bool = True
    RhieChowScaleFactor: float | int = 1.0
    multiMaterialCorrection: multi_material.MultiMaterialInterface = field(factory=lambda: multi_material.Uniform(defaultWeights=1, defaultWeightsGrad=1))

Constant = MechanicsSubSolver  # alias

@ffn_define
class LargeStrainTotLag(Constant):
    """
    Mechanics solver class that solves for total displacement `D` in a total
    Lagrangian framework (i.e. the mesh is NOT updated), employing a finite-strain
    mechanics approximation.

    Options (inherited)
    -------------------
    forceSummary : bool
        Print a summary of boundary-forces at each iteration.
        (default: False; required: False)

    cylindricalStress : bool
        Creates cylindrical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    sphericalStress : bool
        Creates spherical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    RhieChowCorrection : bool
        Apply a Rhie–Chow–like correction to stabilize the momentum equation.
        (default: True; required: False)

    RhieChowScaleFactor : scalar
        Scaling factor for the Rhie–Chow correction.
        (default: 1; required: False)

    multiMaterialCorrection : dict
        Optional sub-dictionary activating a custom discretization for
        multi-material interfaces.
        (required: False)

    Options
    -------
    strainTensor : word
        Select strain tensor definition.
        (default: 'EulerAlmansi'; required: False)
    """
    TYPE: ClassVar[str] = 'largeStrainTotLag'
    strainTensor: str = 'EulerAlmansi'

@ffn_define
class LargeStrainUpdLag(Constant):
    """
    Mechanics solver class that solves for incremental displacement `DD` in an
    updated Lagrangian configuration (i.e. the mesh is updated at every time step)
    and employs a finite-strain mechanics framework.

    Geometric non-linearity is handled via the relative deformation gradient
    and its Jacobian.

    Options (inherited)
    -------------------
    forceSummary : bool
        Print a summary of boundary-forces at each iteration.
        (default: False; required: False)

    cylindricalStress : bool
        Creates cylindrical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    sphericalStress : bool
        Creates spherical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    RhieChowCorrection : bool
        Apply a Rhie–Chow–like correction to stabilize the momentum equation.
        (default: True; required: False)

    RhieChowScaleFactor : scalar
        Scaling factor for the Rhie–Chow correction.
        (default: 1; required: False)

    multiMaterialCorrection : dict
        Optional sub-dictionary activating a custom discretization for
        multi-material interfaces.
        (required: False)

    Options
    -------
    strainTensor : word
        Select strain tensor definition.
        (default: 'EulerAlmansi'; required: False)
    """
    TYPE: ClassVar[str] = 'largeStrainUpdLag'
    strainTensor: str = 'EulerAlmansi'

@ffn_define
class SmallStrain(Constant):
    """
    Mechanics solver class that solves for total displacement `D` in a total
    Lagrangian framework (i.e. the mesh is NOT updated), and employs the
    small-strain approximation.

    Options (inherited)
    -------------------
    forceSummary : bool
        Print a summary of boundary-forces at each iteration.
        (default: False; required: False)

    cylindricalStress : bool
        Creates cylindrical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    sphericalStress : bool
        Creates spherical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    RhieChowCorrection : bool
        Apply a Rhie–Chow–like correction to stabilize the momentum equation.
        (default: True; required: False)

    RhieChowScaleFactor : scalar
        Scaling factor for the Rhie–Chow correction.
        (default: 1; required: False)

    multiMaterialCorrection : dict
        Optional sub-dictionary activating a custom discretization for
        multi-material interfaces.
        (required: False)
    """
    TYPE: ClassVar[str] = 'smallStrain'

@ffn_define
class SmallStrainIncrementalUpdated(Constant):
    """
    Mechanics solver class that solves for incremental displacement `DD` in
    an updated framework (i.e. the mesh is updated at the start of each time step), while
    employing the small-strain approximation.

    Options (inherited)
    -------------------
    forceSummary : bool
        Print a summary of boundary-forces at each iteration.
        (default: False; required: False)

    cylindricalStress : bool
        Creates cylindrical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    sphericalStress : bool
        Creates spherical stress, strain and displacement fields.
        Useful for visualization and data analysis.
        (default: False; required: False)

    RhieChowCorrection : bool
        Apply a Rhie–Chow–like correction to stabilize the momentum equation.
        (default: True; required: False)

    RhieChowScaleFactor : scalar
        Scaling factor for the Rhie–Chow correction.
        (default: 1; required: False)

    multiMaterialCorrection : dict
        Optional sub-dictionary activating a custom discretization for
        multi-material interfaces.
        (required: False)

    Options
    -------
    updateMesh : bool
        Decide whether the mesh is updated at every time step.
        If false, the solver behaves as a small-strain incremental solver.
        (default: True; required: False)
    """
    TYPE: ClassVar[str] = 'smallStrainIncrementalUpdated'
    updateMesh: bool = True
