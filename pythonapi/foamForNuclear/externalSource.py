from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import addParameter
from foamForNuclear.openfoamFile import OpenFOAMFile
from foamForNuclear.timeProfile import TimeProfile


_EXTERNAL_SOURCE_MODE_TYPES = {"transient", "powerMonitoring"}


class ExternalSource(OpenFOAMFile):
    """
    External source object to simulate for example an ADS.

    You also need to specify `nuFission` and `energyPerFission` in nuclear data.

    Parameters
    ----------
    externalSourceMode : str
        External source mode (options: `transient` or `powerMonitoring`)
    beamEnergy : float
        Energy of the source particle (e.g proton energy for spallation) in
        J/source particle
    nuSource : float
        Source neutron yield per source particle (e.g n/p for spallation
        process)
    externalSourceModulationTimeProfile : TimeProfile
        Time profile representing the source modulation. Used only in
        `transient` mode
    modulationFactor : float
        Only used in `powerMonitoring` mode. Modulation factor (default `0`).
    powerModulationDelay : float
        Only used in `powerMonitoring` mode. Power modulation delay
        (default `0`).
    """

    def __init__(
            self,
            region = "",
            externalSourceMode: str="transient",
            beamEnergy: float=0,
            nuSource: float=0,
            externalSourceModulationTimeProfile: TimeProfile=TimeProfile("table"),
            modulationFactor: float=0,
            powerModulationDelay: float=0
        ):
        super().__init__("externalSource", folder="constant", region=region)

        self.externalSourceMode = externalSourceMode
        self.beamEnergy = beamEnergy
        self.nuSource = nuSource

        # Transient
        self.externalSourceModulationTimeProfile: TimeProfile = externalSourceModulationTimeProfile

        # Power monitoring
        self.modulationFactor = modulationFactor
        self.powerModulationDelay = powerModulationDelay


    @property
    def externalSourceMode(self):
        return self._externalSourceMode

    @externalSourceMode.setter
    def externalSourceMode(self, externalSourceMode):
        check_type("externalSourceMode", externalSourceMode, str)
        if (externalSourceMode is not None):
            check_value("externalSourceMode", externalSourceMode, _EXTERNAL_SOURCE_MODE_TYPES)
        self._externalSourceMode = externalSourceMode


    @property
    def beamEnergy(self):
        return self._beamEnergy

    @beamEnergy.setter
    def beamEnergy(self, beamEnergy):
        check_type("beamEnergy", beamEnergy, (int, float))
        self._beamEnergy = beamEnergy


    @property
    def nuSource(self):
        return self._nuSource

    @nuSource.setter
    def nuSource(self, nuSource):
        check_type("nuSource", nuSource, (int, float))
        self._nuSource = nuSource


    @property
    def modulationFactor(self):
        return self._modulationFactor

    @modulationFactor.setter
    def modulationFactor(self, modulationFactor):
        check_type("modulationFactor", modulationFactor, (int, float))
        self._modulationFactor = modulationFactor


    @property
    def powerModulationDelay(self):
        return self._powerModulationDelay

    @powerModulationDelay.setter
    def powerModulationDelay(self, powerModulationDelay):
        check_type("powerModulationDelay", powerModulationDelay, (int, float))
        self._powerModulationDelay = powerModulationDelay


    @property
    def externalSourceModulationTimeProfile(self):
        return self._externalSourceModulationTimeProfile

    @externalSourceModulationTimeProfile.setter
    def externalSourceModulationTimeProfile(self, externalSourceModulationTimeProfile):
        check_type("externalSourceModulationTimeProfile", externalSourceModulationTimeProfile, TimeProfile)
        self._externalSourceModulationTimeProfile = externalSourceModulationTimeProfile


    @property
    def is_empty(self) -> bool:
        if (self.externalSourceMode == "transient"):
            return(self.externalSourceModulationTimeProfile.is_empty)
        return(False)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""

        text += addParameter('externalSourceMode', self.externalSourceMode, isAddExtraLine=True)

        text += addParameter('beamEnergy', self.beamEnergy, isAddExtraLine=True)
        text += addParameter('nuSource', self.nuSource, isAddExtraLine=True)

        if (self.externalSourceMode == "transient"):
            text += f'externalSourceModulationTimeProfile{self.externalSourceModulationTimeProfile!r}'

        if (self.externalSourceMode == "powerMonitoring"):
            text += addParameter('modulationFactor', self.modulationFactor, isAddExtraLine=True)
            text += addParameter('powerModulationDelay', self.powerModulationDelay, isAddExtraLine=True)

        return(text)
