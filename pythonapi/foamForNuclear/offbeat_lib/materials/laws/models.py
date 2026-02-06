# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict

from . import yield_stress
from . import creep

@ffn_define
class ConstitutiveLaw(FoamForNuclearDict):
    """
    Base constitutive law class.
    """
    TYPE: ClassVar[str] = 'constitutiveLaw'

@ffn_define
class Elasticity(ConstitutiveLaw):
    """
    Class for linear elasticity constitutive law.
    """
    TYPE: ClassVar[str] = 'elasticity'

@ffn_define
class FastNeutronGPLS(ConstitutiveLaw):
    """
    Cconstitutive law for the `fastNeutron`-GPLS model (Guisson–Poussard–Le Saux)
    for Zircaloy-4 under RIA loading conditions (also used in ALCYONE; not intended
    for base irradiation).

    Elastic behaviour is isotropic (Hooke’s law). Viscoplastic flow accounts for
    the orthotropic nature of Zircaloy-4 using the Hill tensor / Hill equivalent
    strain, with no stress threshold between elastic and viscoplastic regimes.

    Irradiation effects are introduced via the fast neutron fluence `phi`, provided
    in units of 1e25 n·m⁻² (i.e. set `phi=1` for 1e25 n·m⁻²). A second parameter `zeta`
    (non-physical) corresponds to the normalised axial position along the full rod
    (0 bottom, 1 top); it weakly affects the strength modulus K (default 0.5).


    Options
    -------
    phi : scalar
        Fast neutron fluence in units of 1e25 n·m⁻² (use 1 for 1e25).
        (required: True)

    zeta : scalar
        Non-physical normalised axial position along the full rod (0 bottom, 1 top).
        (default: 0.5; required: False)

    error : scalar
        Newton–Raphson tolerance.
        (default: '1e-8'; required: False)

    max_iter : label
        Maximum Newton–Raphson iterations.
        (default: 2500; required: False)

    relax : scalar
        Relaxation factor for Newton–Raphson iterations.
        (default: 0.1; required: False)

    theta : scalar
        Integration scheme parameter (0..1):
        0 explicit (Forward Euler), 1 implicit (Backward Euler), 0.5 midpoint / Crank–Nicolson.
        (default: 1; required: False)
    """
    TYPE: ClassVar[str] = 'fastNeutronGPLS'
    phi: float | int
    zeta: float | int = 0.5
    error: float | int = 1e-8
    max_iter: int = 2500
    relax: float | int = 0.1
    theta: float | int = 1.0

@ffn_define
class HydrogenGPLS(ConstitutiveLaw):
    """
    hydrogenGPLS (Guisson–Poussard–Le Saux) constitutive law for Zircaloy-4
    under RIA loading conditions (used in ALCYONE; not intended for base irradiation).

    Elastic behaviour is isotropic (Hooke’s law). Viscoplastic flow accounts for
    the orthotropic nature of irradiated Zircaloy-4 using the Hill tensor and Hill
    equivalent strain, with no stress threshold between elastic and viscoplastic
    regimes.

    Irradiation effects are introduced via the fast neutron fluence `phi`, provided
    in units of 1e25 n·m⁻² (i.e. set `phi=1` for 1e25 n·m⁻²). A second parameter `zeta`
    (non-physical) corresponds to the normalised axial position along the full fuel
    rod (0 bottom, 1 top); it weakly affects the strength modulus K (default 0.5).

    Validity range (as documented):
    - fluence up to 10e25 n·m⁻² (≈ 64 GWd/tU)
    - 20 °C < T < 1100 °C (may be exceeded in RIA simulations)
    - 3e-4 < dε/dt < 5 s⁻¹


    Options
    -------
    phi : scalar
        Fast neutron fluence in units of 1e25 n·m⁻² (use 1 for 1e25).
        (required: True)

    zeta : scalar
        Non-physical normalised axial position along the full rod (0 bottom, 1 top).
        (default: 0.5; required: False)

    error : scalar
        Newton–Raphson tolerance.
        (default: '1e-6'; required: False)

    max_iter : label
        Maximum Newton–Raphson iterations.
        (default: 10000; required: False)

    relax : scalar
        Relaxation factor for Newton–Raphson iterations.
        (default: 0.01; required: False)

    theta : scalar
        Integration scheme parameter (0..1):
        0 explicit (Forward Euler), 1 implicit (Backward Euler), 0.5 midpoint / Crank–Nicolson.
        Fully explicit integration is not advised.
        (default: 1; required: False)
    """
    TYPE: ClassVar[str] = 'hydrogenGPLS'
    phi: float | int
    zeta: float | int = 0.5
    error: float | int = 1e-6
    max_iter: int = 10000
    relax: float | int = 0.01
    theta: float | int = 1.0

@ffn_define
class HyperElasticity(ConstitutiveLaw):
    """
    St. Venant–Kirchhoff **hyperelastic constitutive law**.

    This model describes finite-strain elastic behaviour using the
    St. Venant–Kirchhoff formulation.
    """
    TYPE: ClassVar[str] = 'hyperElasticity'

@ffn_define
class MisesPlasticity(ConstitutiveLaw):
    """
    Linear elastic–Mises plasticity constitutive law.

    Elastic behaviour is linear isotropic. Plastic flow is governed by a
    von Mises yield criterion.

    This law requires a `yieldStress` model to be provided as a sub-dictionary.


    Options
    -------
    yieldStress : yieldStress
        Yield stress model controlling the onset of plastic flow.
        This entry selects a run-time-selectable `yieldStress` class and
        provides its associated parameters.
        (required: True)
    """
    TYPE: ClassVar[str] = 'misesPlasticity'
    yieldStress: yield_stress.YieldStress = field(factory=yield_stress.YieldStress)

@ffn_define
class NeoHookeanElasticity(ConstitutiveLaw):
    """
    Neo-Hookean **elastic constitutive law**.

    This model describes finite-strain elastic behaviour using the
    Neo-Hookean formulation.
    """
    TYPE: ClassVar[str] = 'neoHookeanElasticity'

@ffn_define
class NeoHookeanMisesPlasticity(ConstitutiveLaw):
    """
    Neo-Hookean hyperelastic–Mises plasticity constitutive law.

    Elastic behaviour follows a Neo-Hookean hyperelastic formulation suitable
    for finite strains. Plastic flow is governed by a von Mises yield criterion.

    This law requires a `yieldStress` model to be provided as a sub-dictionary.


    Options
    -------
    yieldStress : yieldStress
        Yield stress model controlling the onset of plastic flow.
        This entry selects a run-time-selectable `yieldStress` class and
        provides its associated parameters.
        (required: True)
    """
    TYPE: ClassVar[str] = 'neoHookeanMisesPlasticity'
    yieldStress: yield_stress.YieldStress = field(factory=yield_stress.YieldStress)

@ffn_define
class HyperElasticMisesPlasticCreep(MisesPlasticity):
    """
    Hyperelastic–Mises plasticity–creep constitutive law.

    Elastic behaviour follows a St. Venant–Kirchhoff hyperelastic formulation,
    suitable for finite strains. Plastic flow is governed by a von Mises yield
    criterion, while time-dependent inelastic deformation is captured through
    an associated creep model.

    This law requires both a `yieldStress` model and a `creep` model to be
    provided as sub-dictionaries.


    Options
    -------
    yieldStress : yieldStress
        Yield stress model controlling the onset of plastic flow.
        This entry selects a run-time-selectable `yieldStress` class and
        provides its associated parameters.
        (required: True)

    creep : creep
        Creep model controlling time-dependent inelastic deformation.
        This entry selects a run-time-selectable `creep` class and provides
        its associated parameters.
        (required: True)
    """
    TYPE: ClassVar[str] = 'hyperElasticMisesPlasticCreep'
    yieldStress: yield_stress.YieldStress = field(factory=yield_stress.YieldStress)
    creep: creep.Creep = field(factory=creep.Creep)

@ffn_define
class MisesPlasticCreep(MisesPlasticity):
    """
    Linear elastic–Mises plasticity–creep constitutive law.

    Elastic behaviour is linear isotropic. Plastic flow is governed by a
    von Mises yield criterion, and time-dependent deformation is captured
    through an associated creep model.


    Options
    -------
    yieldStress : yieldStress
        Yield stress model controlling the onset of plastic flow.
        This entry selects a run-time-selectable `yieldStress` class and
        provides its associated parameters.
        (required: True)

    creep : creep
        Creep model controlling time-dependent inelastic deformation.
        This entry selects a run-time-selectable `creep` class and provides
        its associated parameters.
        (required: True)
    """
    TYPE: ClassVar[str] = 'misesPlasticCreep'
    yieldStress: yield_stress.YieldStress = field(factory=yield_stress.YieldStress)
    creep: creep.Creep = field(factory=creep.Creep)
