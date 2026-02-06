# AUTO-GENERATED FILE. DO NOT EDIT BY HAND.
# Generated from OFFBEAT/OpenFOAM YAML docs + C++ TypeName/inheritance.

from __future__ import annotations
from typing import Any, ClassVar
from attrs import field
from foamForNuclear._attrs_tools import ffn_define
from foamForNuclear.common import FoamForNuclearDict

@ffn_define
class Burnup(FoamForNuclearDict):
    """
    Base burnup model class. It disables burnup tracking entirely, i.e. when this
    model is selected, the burnup field `Bu` is **not created** and is not
    present in the simulation. No burnup evolution or burnup-dependent quantities
    are computed.

    This is the default choice when no burnup model is explicitly selected.

    Other burnup models are available; see the documentation of the corresponding
    burnup implementations for details.
    """
    TYPE: ClassVar[str] = 'none'

@ffn_define
class Constant(Burnup):
    """
    Burnup model that allows the burnup field to be prescribed
    explicitly from a `Bu` field file located in the starting time folder.

    If a `Bu` file is present in the starting time folder, the internal field and
    boundary conditions are read directly from that file. If the `Bu` file is not
    present, the burnup field is initialized to 0 MWd/MT in all cells,
    with `zeroGradient` boundary conditions applied to all non-empty and non-wedge
    patches.

    When this model is selected, the burnup field **does not evolve over time** and
    remains constant throughout the simulation, either as defined in the starting
    time folder or as set by default by OFFBEAT.
    """
    TYPE: ClassVar[str] = 'constant'

@ffn_define
class FromPower(Constant):
    """
    Burnup model that tracks the evolution of the local burnup by integrating the
    the volumetric power density field (divided by local density field) over time.

    The burnup is expressed in MWd/MT$_{oxide}$** and is updated at each time step.


    Options
    -------
    heatSourceName : str
        Name of the volumetric power density field used to compute burnup.
        By default, OFFBEAT expects the power density field to be named `Q`. This
        option might be useful when the power field is provided by an external neutronics
        solver or a custom heat source model using a different field name.
        (default: 'Q'; required: False)
    """
    TYPE: ClassVar[str] = 'fromPower'
    heatSourceName: str = 'Q'

@ffn_define
class Lassmann(FromPower):
    """
    Burnup model that is inspired by Lassmann et al., *“The radial distribution
    of plutonium in high burnup UO2 fuels”* (Journal of Nuclear Materials).

    This model extends the `fromPower` burnup approach by acting as a simplified
    depletion module: it tracks a limited set of relevant nuclides and solves a
    simplified form of the Bateman equations for their time evolution.

    The main purpose is to provide a more accurate radial power distribution in the
    fuel by accounting for spatial variations of nuclide concentrations and neutron
    flux across the radius (per axial slice).

    Options (inherited)
    -------------------
    heatSourceName : str
        Name of the volumetric power density field used to compute burnup.
        By default, OFFBEAT expects the power density field to be named `Q`. This
        option might be useful when the power field is provided by an external neutronics
        solver or a custom heat source model using a different field name.
        (default: 'Q'; required: False)

    Options
    -------
    convergencePrecision : float
        Convergence precision for the depletion iterations in the Lassmann burnup
        calculations. The suggested default is `1e-2`.
        (default: '1e-2'; required: False)
    """
    TYPE: ClassVar[str] = 'Lassmann'
    convergencePrecision: float = 1e-2

@ffn_define
class LassmannFBR(FromPower):
    """
    Burnup model that extends the similar `Lassmann` depletion approach to
    fast-spectrum reactor applications.

    The model is inspired by the work of Lassmann et al.,
    *“The radial distribution of plutonium in high burnup UO2 fuels”*
    (Journal of Nuclear Materials), and generalizes the original formulation to
    mixed-oxide and fast-reactor conditions.

    Compared to the standard `Lassmann` model, this extension includes:

    - Depletion of an extended set of nuclides:
    He4, U234, U235, U236, U237, U238, Np237, Np238, Np239, Pu238, Pu239, Pu240,
    Pu241, Pu242, Pu243, Am241, Am242, Am242m, Am243, Am244, Cm242, Cm243, Cm244,
    Cm245.
    - Explicit accounting for helium production due to:
    alpha decay, ternary fission, and the O-16 (n,α) reaction.
    - Use of parametrized microscopic cross sections read from lookup tables
    located in `constant/XSdata`.

    The parametrization of cross sections depends on the reactor type specified
    by the user:

    - For `LWR`, fission and capture cross sections are parametrized as a function
    of the U-235 enrichment in UO2 fuel.
    - For `FBR`, fission and capture cross sections are parametrized as a function
    of the plutonium weight fraction in MOX fuel.

    Options (inherited)
    -------------------
    heatSourceName : str
        Name of the volumetric power density field used to compute burnup.
        By default, OFFBEAT expects the power density field to be named `Q`. This
        option might be useful when the power field is provided by an external neutronics
        solver or a custom heat source model using a different field name.
        (default: 'Q'; required: False)

    Options
    -------
    convergencePrecision : float
        Convergence precision for the depletion iterations in the Lassmann burnup
        calculations. The suggested default is `1e-2`.
        (default: '1e-2'; required: False)
    """
    TYPE: ClassVar[str] = 'LassmannFBR'
    convergencePrecision: float = 1e-2
