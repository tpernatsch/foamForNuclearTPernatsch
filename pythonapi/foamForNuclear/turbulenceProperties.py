from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.openfoamFile import OpenFOAMFile


_TURBULENCE_TYPE_TYPES = {"laminar", "RAS", "LES"}
_RAS_MODEL_TYPES = {"porousKEpsilon", "laminar", "kEpsilon"}


class RASoptionsObject(OpenFOAMDict):
    """
    Parameters
    ----------
    RASModel : str {"porousKEpsilon", "laminar", "kEpsilon"}
        RAS model name.
    turbulence : bool
        default `True`
    printCoeffs : bool
        default `True`
    """

    def __init__(
            self,
            RASModel: str=None,
            turbulence: bool=True,
            printCoeffs: bool=True
        ):
        super().__init__()
        self.RASModel = RASModel
        self.turbulence = turbulence
        self.printCoeffs = printCoeffs


    def __repr__(self, depth = 0):
        return super().__repr__(depth)


    @property
    def RASModel(self):
        return self._RASModel

    @RASModel.setter
    def RASModel(self, RASModel) -> None:
        check_type("RASModel", RASModel, str, none_ok=True)
        if (RASModel is not None):
            check_value("RASModel", RASModel, _RAS_MODEL_TYPES)
        self._RASModel = RASModel
        self.__setitem__('RASModel', self.RASModel)

    @property
    def turbulence(self):
        return self._turbulence

    @turbulence.setter
    def turbulence(self, turbulence) -> None:
        check_type("turbulence", turbulence, bool)
        self._turbulence = turbulence
        self.__setitem__('turbulence', self.turbulence)

    @property
    def printCoeffs(self):
        return self._printCoeffs

    @printCoeffs.setter
    def printCoeffs(self, printCoeffs) -> None:
        check_type("printCoeffs", printCoeffs, bool)
        self._printCoeffs = printCoeffs
        self.__setitem__('printCoeffs', self.printCoeffs)


class PorousKEpsilonPropertiesPerZone(OpenFOAMDict):
    """
    Parameters
    ----------
    zones : list[str]
        Name of the cellZone to apply the turbluence parameters
    convergenceLength : float
        k and epsilon will exponentially converge to equilibrium according to
        this exponent (default `0.1`)
    turbulenceIntensityCoeff : float
        Turbulent intensity correlation in the form 0.16*Reynolds^-0.125
        Reynolds number calculated by the thermal-hydraulic class, according
        to the input data in phaseProperties (default `0.16`)
    turbulenceIntensityExp : float
        Turbulent intensity correlation in the form 0.16*Reynolds^-0.125
        Reynolds number calculated by the thermal-hydraulic class, according
        to the input data in phaseProperties (default `-0.125`)
    turbulenceLengthScaleCoeff : float
        L = 0.07*Dh (Dh is the hydraulic diameter secificied in phaseProperties)
        (default `0.07`)
    DhStruct : float
        Hydraulic diameter (default `1`)
    """

    def __init__(
            self,
            zones: list[str],
            convergenceLength: float=0.1,
            turbulenceIntensityCoeff: float=0.16,
            turbulenceIntensityExp: float=-0.125,
            turbulenceLengthScaleCoeff: float=0.07,
            DhStruct: float=1,
        ):
        super().__init__()
        self.zones = zones
        self.convergenceLength = convergenceLength
        self.turbulenceIntensityCoeff = turbulenceIntensityCoeff
        self.turbulenceIntensityExp = turbulenceIntensityExp
        self.turbulenceLengthScaleCoeff = turbulenceLengthScaleCoeff
        self.DhStruct = DhStruct


    def __repr__(self, depth = 0):
        textZones = "\"" + ':'.join(self.zones) + "\""

        return textZones + super().__repr__(depth)


    @property
    def zones(self):
        return self._zones

    @zones.setter
    def zones(self, zones) -> None:
        check_type("zones", zones, list, none_ok=True)
        self._zones = zones

    @property
    def convergenceLength(self):
        return self._convergenceLength

    @convergenceLength.setter
    def convergenceLength(self, convergenceLength) -> None:
        check_type("convergenceLength", convergenceLength, (float, int))
        self._convergenceLength = convergenceLength
        self.__setitem__('convergenceLength', self.convergenceLength)

    @property
    def turbulenceIntensityCoeff(self):
        return self._turbulenceIntensityCoeff

    @turbulenceIntensityCoeff.setter
    def turbulenceIntensityCoeff(self, turbulenceIntensityCoeff) -> None:
        check_type("turbulenceIntensityCoeff", turbulenceIntensityCoeff, (float, int))
        self._turbulenceIntensityCoeff = turbulenceIntensityCoeff
        self.__setitem__('turbulenceIntensityCoeff', self.turbulenceIntensityCoeff)

    @property
    def turbulenceIntensityExp(self):
        return self._turbulenceIntensityExp

    @turbulenceIntensityExp.setter
    def turbulenceIntensityExp(self, turbulenceIntensityExp) -> None:
        check_type("turbulenceIntensityExp", turbulenceIntensityExp, (float, int))
        self._turbulenceIntensityExp = turbulenceIntensityExp
        self.__setitem__('turbulenceIntensityExp', self.turbulenceIntensityExp)

    @property
    def turbulenceLengthScaleCoeff(self):
        return self._turbulenceLengthScaleCoeff

    @turbulenceLengthScaleCoeff.setter
    def turbulenceLengthScaleCoeff(self, turbulenceLengthScaleCoeff) -> None:
        check_type("turbulenceLengthScaleCoeff", turbulenceLengthScaleCoeff, (float, int))
        self._turbulenceLengthScaleCoeff = turbulenceLengthScaleCoeff
        self.__setitem__('turbulenceLengthScaleCoeff', self.turbulenceLengthScaleCoeff)

    @property
    def DhStruct(self):
        return self._DhStruct

    @DhStruct.setter
    def DhStruct(self, DhStruct) -> None:
        check_type("DhStruct", DhStruct, (float, int))
        self._DhStruct = DhStruct
        self.__setitem__('DhStruct', self.DhStruct)



class TurbulenceProperties(OpenFOAMFile):
    """
    Turbulence properties object

    Parameters
    ----------
    simulationType : str
        Type of turbulence model (`laminar`, `RAS`, `LES`)
    """

    def __init__(
            self,
            simulationType: str='laminar',
            RASoptions: RASoptionsObject=RASoptionsObject(),
            porousKEpsilonProperties: OpenFOAMListDict=OpenFOAMListDict(PorousKEpsilonPropertiesPerZone, "porousKEpsilonProperties"),
            region: str="",
            ext: str=""
        ):
        super().__init__("turbulenceProperties", "constant", region=region, ext=ext)

        self.simulationType = simulationType
        self.RASoptions = RASoptions
        self.porousKEpsilonProperties = porousKEpsilonProperties


    @property
    def simulationType(self):
        return self._simulationType

    @simulationType.setter
    def simulationType(self, simulationType) -> None:
        check_type("simulationType", simulationType, str)
        check_value("simulationType", simulationType, _TURBULENCE_TYPE_TYPES)
        self._simulationType = simulationType

    @property
    def RASoptions(self):
        return self._RASoptions

    @RASoptions.setter
    def RASoptions(self, RASoptions) -> None:
        check_type("RASoptions", RASoptions, RASoptionsObject)
        self._RASoptions = RASoptions

    @property
    def porousKEpsilonProperties(self):
        return self._porousKEpsilonProperties

    @porousKEpsilonProperties.setter
    def porousKEpsilonProperties(self, porousKEpsilonProperties) -> None:
        check_type("porousKEpsilonProperties", porousKEpsilonProperties, OpenFOAMListDict)
        self._porousKEpsilonProperties = porousKEpsilonProperties


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""
        text += addParameter("simulationType", self.simulationType, isAddExtraLine=True)

        if (self.simulationType == "RAS"):
            text += f"RAS{self.RASoptions!r}\n"

        if (self.simulationType == "RAS" and self.RASoptions.RASModel == "porousKEpsilon"):
            text += f"{self.porousKEpsilonProperties!r}\n"

        return(text)
