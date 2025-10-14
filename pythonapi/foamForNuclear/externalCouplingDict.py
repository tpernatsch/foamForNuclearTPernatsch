from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict, OpenFOAMListDict, Vector, addParameter
from foamForNuclear.openfoamFile import OpenFOAMFile


_EXTERNAL_FUNCTION_TYPES = {
    "fieldIntegralToFMU", "massFlowToFMU", "extSensor", "extPatch"
}


class ExternalFunction(OpenFOAMDict):
    """
    Base class to define External Function objects used for FMU simulations.

    Parameters
    ----------
    name : str
        Name of the ExternalFunction
    type : str
        Type of the ExternalFunction object (`fieldIntegralToFMU`, `extSensor`,
        `extPatch`)
    region : str
        Name of the region
    """

    def __init__(
            self,
            name: str,
            type: str,
            region: str=None,
        ):
        super().__init__(name=name)
        self.type = type
        self.region = region


    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name) -> None:
        check_type("name", name, str)
        self._name = name

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _EXTERNAL_FUNCTION_TYPES)
        self._type_ = type
        self.__setitem__("type", type)

    @property
    def region(self):
        return self._region

    @region.setter
    def region(self, region) -> None:
        check_type("region", region, str, none_ok=True)
        self._region = region
        if (region is not None):
            self.__setitem__("region", region)


class FieldIntegralToFMU(ExternalFunction):
    """
    Output an integral value over a cellZone to an FMU

    Parameters
    ----------
    name : str
        Name of the external function
    nameFMU : str
        Name of the output scalar
    fieldName : str
        Name of the measurement field
    cellZone : str
        Name of the cellZone
    region : str
        Name of the region (default `None`).
    scaleFactor: float
        Scaling factor (default `1`).
    initValue : float
        Initial value (default `0`).
    """
    def __init__(
            self,
            name: str,
            nameFMU: str,
            fieldName: str,
            cellZone: str,
            region: str=None,
            scaleFactor: float=1,
            initValue: float=0
        ):
        super().__init__(name=name, type="fieldIntegralToFMU", region=region)

        self.nameFMU = nameFMU
        self.fieldName = fieldName
        self.cellZone = cellZone
        self.scaleFactor = scaleFactor
        self.initValue = initValue

    @property
    def nameFMU(self):
        return self._nameFMU

    @nameFMU.setter
    def nameFMU(self, nameFMU) -> None:
        check_type("nameFMU", nameFMU, str)
        self._nameFMU = nameFMU
        self.__setitem__("nameFMU", nameFMU)

    @property
    def fieldName(self):
        return self._fieldName

    @fieldName.setter
    def fieldName(self, fieldName) -> None:
        check_type("fieldName", fieldName, str)
        self._fieldName = fieldName
        self.__setitem__("fieldName", fieldName)

    @property
    def cellZone(self):
        return self._cellZone

    @cellZone.setter
    def cellZone(self, cellZone) -> None:
        check_type("cellZone", cellZone, str)
        self._cellZone = cellZone
        self.__setitem__("cellZone", cellZone)

    @property
    def scaleFactor(self):
        return self._scaleFactor

    @scaleFactor.setter
    def scaleFactor(self, scaleFactor) -> None:
        check_type("scaleFactor", scaleFactor, (float, int))
        self._scaleFactor = scaleFactor
        self.__setitem__("scaleFactor", scaleFactor)

    @property
    def initValue(self):
        return self._initValue

    @initValue.setter
    def initValue(self, initValue) -> None:
        check_type("initValue", initValue, (float, int))
        self._initValue = initValue
        self.__setitem__("initValue", initValue)


_REGION_TYPES_TYPES = {"patch", "faceSet", "faceZone"}


class MassFlowToFMU(ExternalFunction):
    """
    Output the mass flow rate at a patch to an FMU

    Parameters
    ----------
    name : str
        Name of the external function
    nameFMU : str
        Name of the output scalar
    region : str
        Name of the region
    regionName : str
        Name of the patch
    regionType : str
        Region type. Options: `"patch"`, `"faceSet"`, `"faceZone"`).
    alphaRhoPhiName : str
        Name of mass face flow field (default `"alphaRhoPhi"`).
    scaleFactor: float
        Scaling factor (default `1`).
    initValue : float
        Initial value (default `0`).
    """
    def __init__(
            self,
            name: str,
            nameFMU: str,
            regionName: str,
            regionType: str,
            alphaRhoPhiName: str="alphaRhoPhi",
            region: str=None,
            scaleFactor: float=1,
            initValue: float=0,
        ):
        super().__init__(name=name, type="massFlowToFMU", region=region)

        self.nameFMU = nameFMU
        self.regionName = regionName
        self.regionType = regionType
        self.alphaRhoPhiName = alphaRhoPhiName
        self.scaleFactor = scaleFactor
        self.initValue = initValue

    @property
    def nameFMU(self):
        return self._nameFMU

    @nameFMU.setter
    def nameFMU(self, nameFMU) -> None:
        check_type("nameFMU", nameFMU, str)
        self._nameFMU = nameFMU
        self.__setitem__("nameFMU", nameFMU)

    @property
    def regionName(self):
        return self._regionName

    @regionName.setter
    def regionName(self, regionName) -> None:
        check_type("regionName", regionName, str)
        self._regionName = regionName
        self.__setitem__("regionName", regionName)

    @property
    def regionType(self):
        return self._regionType

    @regionType.setter
    def regionType(self, regionType) -> None:
        check_type("regionType", regionType, str)
        check_value("regionType", regionType, _REGION_TYPES_TYPES)
        self._regionType = regionType
        self.__setitem__("regionType", regionType)

    @property
    def alphaRhoPhiName(self):
        return self._alphaRhoPhiName

    @alphaRhoPhiName.setter
    def alphaRhoPhiName(self, alphaRhoPhiName) -> None:
        check_type("alphaRhoPhiName", alphaRhoPhiName, str)
        self._alphaRhoPhiName = alphaRhoPhiName
        self.__setitem__("alphaRhoPhiName", alphaRhoPhiName)

    @property
    def scaleFactor(self):
        return self._scaleFactor

    @scaleFactor.setter
    def scaleFactor(self, scaleFactor) -> None:
        check_type("scaleFactor", scaleFactor, (float, int))
        self._scaleFactor = scaleFactor
        self.__setitem__("scaleFactor", scaleFactor)

    @property
    def initValue(self):
        return self._initValue

    @initValue.setter
    def initValue(self, initValue) -> None:
        check_type("initValue", initValue, (float, int))
        self._initValue = initValue
        self.__setitem__("initValue", initValue)


class ExtSensor(ExternalFunction):
    """
    FMI output as a single value probing a point in the mesh

    Parameters
    ----------
    name : str
        Name of the external function
    fieldName : str
        Name of the measurement field
    sensorName : str
        Name of the sensor (FMI port name)
    sensorPosition : Vector
        Position of the sensor
    region : str
        Name of the region
    initValue : float
        Initial value (default `None`).
    """

    def __init__(
            self,
            name: str,
            fieldName: str,
            sensorName: str,
            sensorPosition: Vector,
            region: str=None,
            initValue: float=None
        ):
        super().__init__(name=name, type="extSensor", region=region)

        self.fieldName = fieldName
        self.sensorName = sensorName
        self.sensorPosition = sensorPosition
        self.initValue = initValue

    @property
    def fieldName(self):
        return self._fieldName

    @fieldName.setter
    def fieldName(self, fieldName) -> None:
        check_type("fieldName", fieldName, str)
        self._fieldName = fieldName
        self.__setitem__("fieldName", fieldName)

    @property
    def sensorName(self):
        return self._sensorName

    @sensorName.setter
    def sensorName(self, sensorName) -> None:
        check_type("sensorName", sensorName, str)
        self._sensorName = sensorName
        self.__setitem__("sensorName", sensorName)

    @property
    def sensorPosition(self):
        return self._sensorPosition

    @sensorPosition.setter
    def sensorPosition(self, sensorPosition) -> None:
        check_type("sensorPosition", sensorPosition, Vector)
        self._sensorPosition = sensorPosition
        self.__setitem__("sensorPosition", sensorPosition)

    @property
    def initValue(self):
        return self._initValue

    @initValue.setter
    def initValue(self, initValue) -> None:
        check_type("initValue", initValue, (int, float), none_ok=True)
        self._initValue = initValue
        if (initValue is not None):
            self.__setitem__("initValue", initValue)


class ExtPatch(ExternalFunction):
    """
    FMI output on a patch

    Parameters
    ----------
    name : str
        Name of the external function
    fieldName : str
        Name of the measurement field
    outputName : str
        Name of the sensor (FMI port name)
    patchName : str
        Patch name
    region : str
        Name of the region
    initValue : float
        Initial value (default `None`).
    """

    def __init__(
            self,
            name: str,
            fieldName: str,
            outputName: str,
            patchName: str,
            region: str=None,
            initValue: float=None
        ):
        super().__init__(name=name, type="extPatch", region=region)

        self.fieldName = fieldName
        self.outputName = outputName
        self.patchName = patchName
        self.initValue = initValue

    @property
    def fieldName(self):
        return self._fieldName

    @fieldName.setter
    def fieldName(self, fieldName) -> None:
        check_type("fieldName", fieldName, str)
        self._fieldName = fieldName
        self.__setitem__("fieldName", fieldName)

    @property
    def outputName(self):
        return self._outputName

    @outputName.setter
    def outputName(self, outputName) -> None:
        check_type("outputName", outputName, str)
        self._outputName = outputName
        self.__setitem__("outputName", outputName)

    @property
    def patchName(self):
        return self._patchName

    @patchName.setter
    def patchName(self, patchName) -> None:
        check_type("patchName", patchName, str)
        self._patchName = patchName
        self.__setitem__("patchName", patchName)

    @property
    def initValue(self):
        return self._initValue

    @initValue.setter
    def initValue(self, initValue) -> None:
        check_type("initValue", initValue, (int, float), none_ok=True)
        self._initValue = initValue
        if (initValue is not None):
            self.__setitem__("initValue", initValue)



class ExternalCouplingDict(OpenFOAMFile):
    """
    Main file containing all the external coupling functions for FMI coupling.

    Parameters
    ----------
    externalFunctions : OpenFOAMListDict
        Default `None`.
    """
    def __init__(
            self,
            externalFunctions: OpenFOAMListDict=None
        ):
        super().__init__("externalCouplingDict", folder="system")

        self.externalFunctions = externalFunctions


    @property
    def externalFunctions(self):
        return self._externalFunctions

    @externalFunctions.setter
    def externalFunctions(self, externalFunctions) -> None:
        if externalFunctions is not None:
            check_type("externalFunctions", externalFunctions, OpenFOAMListDict)
            self._externalFunctions = externalFunctions
        else:
            self._externalFunctions = OpenFOAMListDict(ExternalFunction, "externalFunctions")

    def append(self, externalFunction: ExternalFunction):
        check_type("externalFunction", externalFunction, ExternalFunction)
        self.externalFunctions.append(externalFunction)


    @OpenFOAMFile._write_to_file
    def export_to_openfoam(self):
        text = ""
        text += f"{self.externalFunctions}\n"

        return(text)
