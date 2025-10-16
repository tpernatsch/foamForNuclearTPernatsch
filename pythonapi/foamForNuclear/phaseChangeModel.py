from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import List, OpenFOAMDict


_PHASE_CHANGE_MODEL_TYPES = {"heatDriven"}
_HEAT_DRIVEN_MODES = {
    "conductionLimited", "onePhaseDriven", "twoPhaseDriven", "mixedDriven"
}
_LATENT_HEAT_TYPES = {"FinkLeibowitz", "fromThermophysicalProperties", "water"}
_SATURATION_MODEL_TYPES = {
    "BrowningPotter", "constantTemperature", "water", "waterTRACE"
}


class LatentHeatModel(OpenFOAMDict):
    def __init__(
            self,
            type: str,
            adjust: bool=False
        ):
        super().__init__()

        self.type = type
        self.adjust = adjust

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _LATENT_HEAT_TYPES)
        self._type_ = type
        self.__setitem__("type", type)

    @property
    def adjust(self):
        return self._adjust

    @adjust.setter
    def adjust(self, adjust) -> None:
        check_type("adjust", adjust, bool)
        self._adjust = adjust
        self.__setitem__("adjust", adjust)


class FinkLeibowitzLatentHeat(LatentHeatModel):
    """
    This class describes the latent heat of vaporization of liquid sodium
    according to the data presented in the report by Fink & Leibowitz
    (ANL RE 95/2), https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf
    """
    def __init__(self, adjust = False):
        super().__init__(type="FinkLeibowitz", adjust=adjust)


class FromThermophysicalPropertiesLatentHeat(LatentHeatModel):
    """
    Class that computes latent heat based on the specified enthalpies of
    formation of the fluids. These are specified under the Hf keyword
    in the thermophysical properties of the two fluids. The latent heat is thus
    computed as Hf.vapour - Hf.liquid.
    """
    def __init__(self, adjust = False):
        super().__init__(type="fromThermophysicalProperties", adjust=adjust)


class WaterLatentHeat(LatentHeatModel):
    """
    This class describes the latent heat of vaporization of water in the
    0.01-350 C range. It was obtained from a fit (in the form of
    L = A+B*log(C-T) with T being the liquid temperature (K)) of NIST water
    data at saturation found at:
    https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf

    It performs reasonably well up to ~ 365 C as well, considering that the
    critical water temperature is 375.16 C (at which point the latent heat
    becomes 0, but is limited for obvious numerical reasons as this code is not
    meant for sub-supercritical fluid transition simulations)
    """
    def __init__(self, adjust = False):
        super().__init__(type="water", adjust=adjust)


class SaturationModel(OpenFOAMDict):
    def __init__(
            self,
            type: str,
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _SATURATION_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class BrowningPotterSaturationModel(SaturationModel):
    """
    Saturation model based on results by Browning and Potter, refer to
    https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf. The saturation pressure form is
    lnPSat = A + B/TSat + C*ln(TSat), which is not invertible in TSat. TSat is
    inverted by approximateding ln(TSat) by a second order polyinomial.
    """
    def __init__(self):
        super().__init__(type="BrowningPotter")


class ConstantTemperature(SaturationModel):
    """
    Constant saturation pressure and temperature.

    Parameters
    ----------
    value : float
        Saturation temperature
    """
    def __init__(self, value: float):
        super().__init__(type="constantTemperature")

        self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        check_type("value", value, (float, int))
        self._value = value
        self.__setitem__("value", value)


class WaterSaturationModel(SaturationModel):
    """
    Saturation model for water in the 0.01-350 C range based on an
    interpolation of NIST data in the form

    .. math::
        Psat = A (T-B)^C

    The NIST data can be found at this `link
    <https://www.nist.gov/system/files/documents/srd/NISTIR5078-Tab1.pdf>`_.
    """
    def __init__(self):
        super().__init__(type="water")


class WaterTRACESaturationModel(SaturationModel):
    """
    Saturation model for water in the 0.01-350 C range based on `TRACE
    <https://www.nrc.gov/docs/ML1200/ML120060218.pdf>`_.
    """
    def __init__(self):
        super().__init__(type="waterTRACE")


class PhaseChangeModel(OpenFOAMDict):
    def __init__(
            self,
            type: str,
            latentHeatModel: LatentHeatModel=None,
            saturationModel: SaturationModel=None,
            correctLatentHeat: bool=False,
            residualInterfacialArea: float=1e-3
        ):
        super().__init__(name="phaseChangeModel")

        self.type = type
        self.correctLatentHeat = correctLatentHeat
        self.residualInterfacialArea = residualInterfacialArea
        self.latentHeatModel = latentHeatModel
        self.saturationModel = saturationModel


    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _PHASE_CHANGE_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)

    @property
    def latentHeatModel(self):
        return self._latentHeatModel

    @latentHeatModel.setter
    def latentHeatModel(self, latentHeatModel) -> None:
        check_type("latentHeatModel", latentHeatModel, LatentHeatModel, none_ok=True)
        self._latentHeatModel = latentHeatModel
        if (latentHeatModel is not None):
            self.__setitem__("latentHeatModel", latentHeatModel)

    @property
    def saturationModel(self):
        return self._saturationModel

    @saturationModel.setter
    def saturationModel(self, saturationModel) -> None:
        check_type("saturationModel", saturationModel, SaturationModel, none_ok=True)
        self._saturationModel = saturationModel
        if (saturationModel is not None):
            self.__setitem__("saturationModel", saturationModel)


class ForcedConstant(PhaseChangeModel):
    """
    """
    def __init__(
            self,
            value: float,
            regions: list[str],
            latentHeatModel: LatentHeatModel=None,
            saturationModel: SaturationModel=None,
            correctLatentHeat: bool=False
        ):
        super().__init__(
            "forcedConstant",
            latentHeatModel=latentHeatModel,
            saturationModel=saturationModel,
            correctLatentHeat=correctLatentHeat
        )

        self.value = value
        self.regions = regions


    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        check_type("value", value, (float, int))
        self._value = value
        self.__setitem__("value", value)

    @property
    def regions(self):
        return self._regions

    @regions.setter
    def regions(self, regions) -> None:
        check_type("regions", regions, list)
        self._regions = regions
        self.__setitem__("regions", List(regions))


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
    def __init__(
            self,
            mode: str,
            drivingPhase: str=None,
            latentHeatModel: LatentHeatModel=None,
            saturationModel: SaturationModel=None,
            correctLatentHeat: bool=False,
            residualInterfacialArea: float=1e-3
        ):
        super().__init__(
            "heatDriven",
            latentHeatModel=latentHeatModel,
            saturationModel=saturationModel,
            correctLatentHeat=correctLatentHeat,
            residualInterfacialArea=residualInterfacialArea
        )

        self.mode = mode
        self.drivingPhase = drivingPhase


    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, mode) -> None:
        check_type("mode", mode, str)
        check_value("mode", mode, _HEAT_DRIVEN_MODES)
        self._mode = mode
        self.__setitem__("mode", mode)

    @property
    def drivingPhase(self):
        return self._drivingPhase

    @drivingPhase.setter
    def drivingPhase(self, drivingPhase) -> None:
        check_type("drivingPhase", drivingPhase, str, none_ok=True)
        self._drivingPhase = drivingPhase
        if (drivingPhase is not None):
            self.__setitem__("drivingPhase", drivingPhase)
