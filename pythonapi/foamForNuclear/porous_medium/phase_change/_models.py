import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OffbeatDict, List, OpenFOAMDict
from foamForNuclear._attrs_tools import offbeat_define, _to_List_str
from typing import ClassVar
from attrs import field, validators as v

from .latent_heat import LatentHeatModel
from .saturation import SaturationModel


@offbeat_define
class PhaseChangeModel(OffbeatDict):
    TYPE: ClassVar[str] = "none"
    latentHeatModel: LatentHeatModel | None = None,
    saturationModel: SaturationModel | None = None,
    correctLatentHeat: bool = False,
    residualInterfacialArea: int | float = 1e-3


@offbeat_define
class ForcedConstant(PhaseChangeModel):
    """
    """
    TYPE: ClassVar[str] = "forcedConstant"
    value: int | float
    regions: list[str] | List[str] = field(factory=list)


_HEAT_DRIVEN_MODES = {"conductionLimited", "twoPhaseDriven", "onePhaseDriven"}

@offbeat_define
class HeatDriven(PhaseChangeModel):
    """
    Phase change models that compute mass transfer based on interfacial
    heat fluxes. There are three possible modes that can be selected:

    -   `conductionLimited`: `dmdt` is computed so to conserve total energy
        transfer across interface if the interfacial heat fluxes do not
        balance out. It q1i and q2i are the respective heat fluxes form the
        bulk of phase 1/2 to the interface, then dmdt = (q1i+q2i)/L with L
        being the latent heat. This is the approach found in TRACE
        (check theory manual, https://www.nrc.gov/docs/ML0710/ML071000097.pdf)
        but it is not suited if there is a massive difference in the volumetric
        heat capacity (J/K/m3) of the phases (such as in the case of sodium, due
        to the density difference of a factor ~ 2000). It can be very unstable
        and not converge at all in such scenarios;

    -   `twoPhaseDriven`: evaporation is driven uniquely by liquid superheat and
        condensation is driven uniquely by vapor undercooling. There is no need
        to specify which phase is liquid and which is vapour in the
        phaseProperties dict, as this information is contained in the latent
        heat L (positive if phase1 is liquid and 2 is vapour, negative
        otherwise);

    -   `onePhaseDriven`: both evaporation and condensation are driven uniquely
        by one phase's superheat or cooling. Implemented as a counterpart to
        twoPhaseDriven.

    The latent heat is computed as a difference of the enthalpies of formation
    (specified in the thermoPhysicalProperties of each phase under the Hf
    keyword). Thus, the vapour phase is the phase with the highest enthalpy of
    formation. The latent heat can be adjusted as described in the TRACE theory
    manual to avoid thermal-run-aways in particular scenarios. This feature can
    be enabled via the correctLatentHeat flag in the phaseChangeModel subDict,
    but it is not recommended if dealing with phases with a low volumetric
    heat capacity, as described before (e.g. Sodium).

    Parameters
    ----------
    mode : str {"conductionLimited", "onePhaseDriven", "twoPhaseDriven", "mixedDriven"}
        Heat driven mode.
    """
    TYPE: ClassVar[str] = "heatDriven"
    mode: str | None = field(
        default = None,
        validator=v.optional(v.and_(v.instance_of(str), v.in_(_HEAT_DRIVEN_MODES),)),
    )
    drivingPhase: str | None = None