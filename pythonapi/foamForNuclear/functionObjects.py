import pandas as pd
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *

# Duplicated from controlDict.py
_WRITE_CONTROL_TYPES = {
    'none', 'timeStep', 'runTime', 'adjustable', 'adjustableRunTime',
    'writeTime', 'clockTime', 'cpuTime'
}


_FUNCTION_OBJECT_TYPES = {
    'fieldMinMax', 'surfaceFieldValue', 'volFieldValue',
    'mag',
    'probes',
    'multiply', 'subtract', 'divide', 'limitFields',
    'massFlow', 'TBulk',
    'FMUSimulator'
}
_FUNCTION_OBJECT_LIBS = {
    "fieldFunctionObjects", "libfieldFunctionObjects.so", "libsampling.so",
    'pyFMUSim', "libFunctionObjects.so", "libFFNFunctionObjects.so",
}
_VOLUME_OPERATION_TYPES = {"volIntegrate"}
_SURFACE_OPERATION_TYPES = {
    "none", "min", "max", "sum", "sumMag", "sumDirection", "sumDirectionBalance",
    "average", "areaAverage", "areaIntegrate", "CoV", "areaNormalAverage",
    "areaNormalIntegrate", "uniformity", "weightedSum", "weightedAverage",
    "weightedAreaAverage", "weightedAreaIntegrate", "weightedUniformity",
    "absWeightedSum", "absWeightedAverage", "absWeightedAreaAverage",
    "absWeightedAreaIntegrate", "absWeightedUniformity"
}
_SURFACE_REGION_TYPES = {
    "faceZone", "patch", "functionObjectSurface", "sampledSurface"
}
_FMI_CONNECTION_DIRECTION_TYPES = {'from', 'to'}

_MAP_METHODS = {"direct", "mapNearest", "cellVolumeWeight", "correctedCellVolumeWeight"}



class FunctionObject(OpenFOAMDict):
    def __init__(
            self,
            name: str,
            type: str,
            libs: str=None,
            log: bool=None,
            writeFields: bool=None,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
            regionType: str=None,
            regionName: str=None,
            alphaRhoPhiName: str=None,
            scaleFactor: float=None
        ):
        super().__init__(name=name)
        self.type = type
        self.libs = libs
        self.log = log
        self.writeFields = writeFields
        self.writeControl = writeControl
        self.writeInterval = writeInterval
        self.region = region
        self.regionType = regionType
        self.regionName = regionName
        self.alphaRhoPhiName = alphaRhoPhiName
        self.scaleFactor = scaleFactor


    def __repr__(self, depth = 0):
        self.__setitem__('libs', f'("{self.libs}")')

        return super().__repr__(depth)


    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type):
        check_type("type", type, str)
        check_value("type", type, _FUNCTION_OBJECT_TYPES)
        self._type_ = type
        if (type in ["surfaceFieldValue"]):
            self.libs = "fieldFunctionObjects"
        self.__setitem__('type', self.type)

    @property
    def libs(self):
        return self._libs

    @libs.setter
    def libs(self, libs):
        if (libs is not None):
            check_type("libs", libs, str)
            check_value("libs", libs, _FUNCTION_OBJECT_LIBS)
            self._libs = libs
        elif (self.libs is None):
            self._libs = None

    @property
    def log(self):
        return self._log

    @log.setter
    def log(self, log):
        check_type("log", log, bool, none_ok=True)
        self._log = log
        if (log is not None):
            self.__setitem__("log", log)

    @property
    def writeFields(self):
        return self._writeFields

    @writeFields.setter
    def writeFields(self, writeFields):
        check_type("writeFields", writeFields, bool, none_ok=True)
        self._writeFields = writeFields
        if (writeFields is not None):
            self.__setitem__("writeFields", writeFields)

    @property
    def writeControl(self):
        return self._writeControl

    @writeControl.setter
    def writeControl(self, writeControl):
        check_type("writeControl", writeControl, str, none_ok=True)
        self._writeControl = writeControl
        if (writeControl is not None):
            check_value("writeControl", writeControl, _WRITE_CONTROL_TYPES)
            self.__setitem__("writeControl", writeControl)

    @property
    def writeInterval(self):
        return self._writeInterval

    @writeInterval.setter
    def writeInterval(self, writeInterval):
        check_type("writeInterval", writeInterval, (float, int), none_ok=True)
        self._writeInterval = writeInterval
        if (writeInterval is not None):
            self.__setitem__("writeInterval", writeInterval)

    @property
    def region(self):
        return self._region

    @region.setter
    def region(self, region):
        check_type("region", region, str, none_ok=True)
        self._region = region
        if (region is not None):
            self.__setitem__("region", region)

    @property
    def regionType(self):
        return self._regionType

    @regionType.setter
    def regionType(self, regionType):
        check_type("regionType", regionType, str, none_ok=True)
        self._regionType = regionType
        if (regionType is not None):
            self.__setitem__("regionType", regionType)

    @property
    def regionName(self):
        return self._regionName

    @regionName.setter
    def regionName(self, regionName):
        check_type("regionName", regionName, str, none_ok=True)
        self._regionName = regionName
        if (regionName is not None):
            self.__setitem__("regionName", regionName)

    @property
    def alphaRhoPhiName(self):
        return self._alphaRhoPhiName

    @alphaRhoPhiName.setter
    def alphaRhoPhiName(self, alphaRhoPhiName):
        check_type("alphaRhoPhiName", alphaRhoPhiName, str, none_ok=True)
        self._alphaRhoPhiName = alphaRhoPhiName
        if (alphaRhoPhiName is not None):
            self.__setitem__("alphaRhoPhiName", alphaRhoPhiName)

    @property
    def scaleFactor(self):
        return self._scaleFactor

    @scaleFactor.setter
    def scaleFactor(self, scaleFactor):
        check_type("scaleFactor", scaleFactor, (float, int), none_ok=True)
        self._scaleFactor = scaleFactor
        if (scaleFactor is not None):
            self.__setitem__("scaleFactor", scaleFactor)

    def read_from_case(self, startTime: float, caseFolder: str='./'):
        msg = "FunctionObject.read_from_case not implemented"
        raise NotImplementedError(msg)


class MapFields(FunctionObject):
    def __init__(
            self,
            name: str,
            fields: list[str],
            mapRegion: str,
            mapMethod: str,
            consitent: bool=True,
            log: bool=None,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
        ):
        """
        Parameters
        ----------
        fields : list[str]
            List of fields name to map from `region` to `mapRegion`.
        mapRegion : str
            Name of the region to map to.
        mapMethod : str {direct, mapNearest, cellVolumeWeight, correctedCellVolumeWeight}
            Mapping method.
        region : str
            Name of the region to map from.
        """
        super().__init__(
            name,
            "mapFields",
            "libfieldFunctionObjects.so",
            log,
            None,
            writeControl,
            writeInterval,
            region,
            None,
            None,
            None,
            1
        )
        self.fields = fields
        self.mapRegion = mapRegion
        self.mapMethod = mapMethod
        self.consitent = consitent

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fields):
        check_type("fields", fields, list)
        self._fields = fields
        self.__setitem__("fields", List(fields))

    @property
    def mapRegion(self):
        return self._mapRegion

    @mapRegion.setter
    def mapRegion(self, mapRegion):
        check_type("mapRegion", mapRegion, str)
        self._mapRegion = mapRegion
        self.__setitem__("mapRegion", mapRegion)

    @property
    def mapMethod(self):
        return self._mapMethod

    @mapMethod.setter
    def mapMethod(self, mapMethod):
        check_type("mapMethod", mapMethod, str)
        check_value("mapMethod", mapMethod, _MAP_METHODS)
        self._mapMethod = mapMethod
        self.__setitem__("mapMethod", mapMethod)

    @property
    def consitent(self):
        return self._consitent

    @consitent.setter
    def consitent(self, consitent):
        check_type("consitent", consitent, bool)
        self._consitent = consitent
        self.__setitem__("consitent", consitent)


class SurfaceFieldValue(FunctionObject):
    def __init__(
            self,
            name: str,
            fields: list[str],
            operation: str,
            log: bool=None,
            writeFields: bool=None,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
            regionType: str=None,
            regionName: str=None,
            scaleFactor: float=1
        ):
        super().__init__(
            name,
            "surfaceFieldValue",
            "libfieldFunctionObjects.so",
            log,
            writeFields,
            writeControl,
            writeInterval,
            region,
            regionType,
            regionName,
            None,
            scaleFactor
        )
        self.fields = fields
        self.operation = operation

    @property
    def regionName(self):
        return self._regionName

    @regionName.setter
    def regionName(self, regionName):
        check_type("regionName", regionName, str)
        self._regionName = regionName
        self.__setitem__("name", regionName)

    @property
    def regionType(self):
        return self._regionType

    @regionType.setter
    def regionType(self, regionType):
        check_type("regionType", regionType, str)
        check_value("regionType", regionType, _SURFACE_REGION_TYPES)
        self._regionType = regionType
        self.__setitem__("regionType", regionType)

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fields):
        check_type("fields", fields, list)
        self._fields = fields
        self.__setitem__("fields", List(fields))

    @property
    def operation(self):
        return self._operation

    @operation.setter
    def operation(self, operation):
        check_type("operation", operation, str)
        check_value("operation", operation, _SURFACE_OPERATION_TYPES)
        self._operation = operation
        self.__setitem__("operation", operation)

    def read_from_case(
            self,
            startTime: float,
            caseFolder: str='./'
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time
        caseFolder : str
            Case folder do read from (default `'./'`).
        """
        filename = f'{caseFolder}/postProcessing/{self.region}/{self.name}/{startTime}/surfaceFieldValue.dat'

        data = pd.read_csv(filename, sep='\t', skiprows=4)
        data.columns = data.columns.str.replace('#','')
        data.columns = data.columns.str.replace(' ','')

        return(data)


class VolFieldValue(FunctionObject):
    def __init__(
            self,
            name,
            fields: list[str],
            operation: str,
            log: bool=None,
            writeFields: bool=None,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
            regionType: str=None,
            regionName: str=None,
            scaleFactor: float=1
        ):
        super().__init__(
            name,
            "volFieldValue",
            "fieldFunctionObjects",
            log,
            writeFields,
            writeControl,
            writeInterval,
            region,
            regionType,
            regionName,
            None,
            scaleFactor
        )
        self.fields = fields
        self.operation = operation

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fields):
        check_type("fields", fields, list)
        self._fields = fields
        self.__setitem__("fields", List(fields))

    @property
    def operation(self):
        return self._operation

    @operation.setter
    def operation(self, operation):
        check_type("operation", operation, str)
        check_value("operation", operation, _VOLUME_OPERATION_TYPES)
        self._operation = operation
        self.__setitem__("operation", operation)

    def read_from_case(
            self,
            startTime: float,
            caseFolder: str='./'
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time
        caseFolder : str
            Case folder do read from (default `'./'`).
        """
        filename = f'{caseFolder}/postProcessing/{self.region}/{self.name}/{startTime}/volFieldValue.dat'

        data = pd.read_csv(filename, sep='\t', skiprows=3)
        data.columns = data.columns.str.replace('#','')
        data.columns = data.columns.str.replace(' ','')

        return(data)


class Probes(FunctionObject):
    """
    Probes functionObject
    """
    def __init__(
            self,
            name,
            fields: list[str],
            enabled: bool,
            probeLocations: list[Vector],
            log: bool=None,
            writeFields: bool=None,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
        ):
        super().__init__(
            name,
            "probes",
            "libsampling.so",
            log,
            writeFields,
            writeControl,
            writeInterval,
            region,
            regionType=None,
            regionName=None,
            alphaRhoPhiName=None,
            scaleFactor=None
        )
        self.fields = fields
        self.enabled = enabled
        self.probeLocations = List([])
        for probeLocation in probeLocations:
            self.add_probe_location(probeLocation)

    def __repr__(self, depth=0):
        self.__setitem__("probeLocations", self.probeLocations)
        return super().__repr__(depth)

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fields):
        check_type("fields", fields, list)
        self._fields = fields
        self.__setitem__("fields", List(fields))

    @property
    def enabled(self):
        return self._enabled

    @enabled.setter
    def enabled(self, enabled):
        check_type("enabled", enabled, bool)
        self._enabled = enabled
        self.__setitem__("enabled", enabled)

    def add_probe_location(self, probeLocation: Vector | tuple) -> None:
        check_type("probeLocation", probeLocation, (Vector, tuple, list, np.ndarray))
        if (isinstance(probeLocation, (tuple, list, np.ndarray))):
            probeLocation = Vector(probeLocation[0], probeLocation[1], probeLocation[2])

        self.probeLocations.append(probeLocation)

    def read_from_case(
            self,
            startTime: float,
            fieldName: str,
            caseFolder: str='./'
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time
        fieldName : str
            Name of the field
        caseFolder : str
            Case folder do read from (default `'./'`).

        Return
        ------
        Return the probe values and the data locations
        """
        region = "" if self.region is None else self.region
        filename = f'{caseFolder}/postProcessing/{self.name}/{region}/{startTime}/{fieldName}'

        locations = []
        data = {}
        with open(filename, 'r') as f:
            for line in f:
                if ("# Probe " in line):
                    pos = line.replace('(', '').replace(')', '').split()[-3:]
                    pos = [float(val) for val in pos]
                    vec = Vector(pos[0], pos[1], pos[2])
                    locations.append(vec)
                elif ("# Time " in line):
                    continue
                else:
                    row = [float(val) for val in line.split()]
                    time, values = row[0], row[1:]
                    data[time] = np.array(values)

        return(data, locations)


class MassFlow(FunctionObject):
    """
    Compute the mass flow rate in kg/s
    """
    def __init__(
            self,
            name,
            log: bool,
            writeFields: bool,
            writeControl: str,
            writeInterval: float,
            region: str,
            regionName: str,
            regionType: str="patch",
            alphaRhoPhiName: str="alphaRhoPhi",
            scaleFactor: float=1
        ):
        super().__init__(
            name,
            type="massFlow",
            libs="libFFNFunctionObjects.so",
            log=log,
            writeFields=writeFields,
            writeControl=writeControl,
            writeInterval=writeInterval,
            region=region,
            regionName=regionName,
            regionType=regionType,
            alphaRhoPhiName=alphaRhoPhiName,
            scaleFactor=scaleFactor
        )

    def read_from_case(
            self,
            startTime: float,
            caseFolder: str='./'
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time
        caseFolder : str
            Case folder do read from (default `'./'`).
        """
        filename = f'{caseFolder}/postProcessing/{self.region}/{self.name}/{startTime}/massFlow.dat'

        data = pd.read_csv(filename, sep=' ', skiprows=2, names=['Time', 'MassFlow'])

        return(data)


class TBulk(FunctionObject):
    """
    Compute the bulk temperature in K.
    """
    def __init__(
            self,
            name,
            log: bool,
            writeFields: bool,
            writeControl: str,
            writeInterval: float,
            region: str,
            regionName: str,
            regionType: str="patch",
            thermoName: str="thermophysicalProperties",
            alphaRhoPhiName: str="alphaRhoPhi",
        ):
        super().__init__(
            name,
            type="TBulk",
            libs="libFFNFunctionObjects.so",
            log=log,
            writeFields=writeFields,
            writeControl=writeControl,
            writeInterval=writeInterval,
            region=region,
            regionName=regionName,
            regionType=regionType,
            alphaRhoPhiName=alphaRhoPhiName
        )
        self.thermoName = thermoName


    @property
    def thermoName(self):
        return self._thermoName

    @thermoName.setter
    def thermoName(self, thermoName):
        check_type("thermoName", thermoName, str, none_ok=True)
        self._thermoName = thermoName
        if (thermoName is not None):
            self.__setitem__("thermoName", thermoName)

    def read_from_case(
            self,
            startTime: float,
            caseFolder: str='./'
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time
        caseFolder : str
            Case folder do read from (default `'./'`).
        """
        filename = f'{caseFolder}/postProcessing/{self.region}/{self.name}/{startTime}/TBulk.dat'

        data = pd.read_csv(filename, sep=' ', skiprows=2, names=['Time', 'TBulk'])

        return(data)


class FMUSimulator(FunctionObject):
    """
    Function object for simulation with FMU.

    Parameters
    ----------
    name : str
        Name of the functionObject
    pyClassName : str
        Name of the Python class to run the FMU simulator. Can be derived from
        `OMSimulatorContainer` or `PyFMIContainer`.
    pyFileName : str
        Name of the Python file to look for the `pyClassName`
    mapping : list[tuple]
        List of tuple with the follow form: `('from'/'to', 'openfoamName', 'fmuName')`.

        Direction of the FMI connection. Options: `'from'` or `'to'`.
        - `'from'` means "FMI value **FROM** the OpenFOAM case"
        - `'to'`   means "FMI value **TO**   the OpenFOAM case"
    """
    def __init__(
            self,
            name,
            pyClassName: str,
            pyFileName: str,
            mapping: list[tuple]
        ):
        super().__init__(name, type="FMUSimulator", libs="pyFMUSim")

        self.pyClassName = pyClassName
        self.pyFileName = pyFileName
        self.mapping = []
        for direction, openfoamName, fmuName in mapping:
            self.add_mapping_connection(direction, openfoamName, fmuName)


    def __repr__(self, depth: int=0):
        self.generate_coupling_file()
        return super().__repr__(depth)

    @property
    def pyClassName(self):
        return self._pyClassName

    @pyClassName.setter
    def pyClassName(self, pyClassName):
        check_type("pyClassName", pyClassName, str)
        self._pyClassName = pyClassName
        self.__setitem__("pyClassName", pyClassName)

    @property
    def pyFileName(self):
        return self._pyFileName

    @pyFileName.setter
    def pyFileName(self, pyFileName):
        check_type("pyFileName", pyFileName, str)
        self._pyFileName = pyFileName
        self.__setitem__("pyFileName", pyFileName)


    def add_mapping_connection(
            self,
            direction: str,
            openfoamName: str,
            fmuName: str
        ) -> None:
        check_type("direction", direction, str)
        check_type("openfoamName", openfoamName, str)
        check_type("fmuName", fmuName, str)
        check_value("direction", direction, _FMI_CONNECTION_DIRECTION_TYPES)

        self.mapping.append((direction, openfoamName, fmuName))


    def get_ports_names(self, direction: str) -> list[str]:
        """
        Parameters
        ----------
        direction : str
            Direction of the FMI connection. Options: `'from'` or `'to'`.
            - `'from'` means "FMI value **FROM** the OpenFOAM case"
            - `'to'`   means "FMI value **TO**   the OpenFOAM case"
        """
        check_type("direction", direction, str)
        check_value("direction", direction, _FMI_CONNECTION_DIRECTION_TYPES)
        return([
            ofName for dir_, ofName, _ in self.mapping if dir_ == direction
        ])


    def generate_coupling_file(self) -> None:
        """
        Generate the FMU.json file used for the mapping between the OpenFOAM
        case and the FMU.
        """
        fromPorts = self.get_ports_names("from")
        toPorts = self.get_ports_names("to")

        text = "{\n"

        text += tab + '"mapping": [\n'
        for i, (_, ofName, fmiName) in enumerate(self.mapping):
            text += f'{2*tab}["{ofName}", "{fmiName}"]{"," if i < len(self.mapping)-1 else ""}\n'
        text += tab + "],\n"

        text += tab + '"from_OF": {\n'
        for i, ofName in enumerate(fromPorts):
            text += f'{2*tab}"{ofName}": "REAL"{"," if i < len(fromPorts)-1 else ""}\n'
        text += tab + "},\n"

        text += tab + '"to_OF": {\n'
        for i, ofName in enumerate(toPorts):
            text += f'{2*tab}"{ofName}": "REAL"{"," if i < len(toPorts)-1 else ""}\n'
        text += tab + "}\n"

        text += "}\n"

        with open("FMU.json", 'w') as f:
            f.write(text)


    def get_coupling_mapping_as_text(self, depth: int=0) -> str:
        text = ""

        text += f"{depth*tab}FMI ports from OpenFOAM:\n"
        for direction, openfoamPort, fmuPort in self.mapping:
            if (direction == "from"):
                text += f"{(depth+1)*tab}{openfoamPort} -> {fmuPort}\n"

        text += f"{depth*tab}FMI ports to OpenFOAM:\n"
        for direction, openfoamPort, fmuPort in self.mapping:
            if (direction == "to"):
                text += f"{(depth+1)*tab}{openfoamPort} <- {fmuPort}\n"

        return(text)


    def plot_coupling_graph(self):
        """
        Generate a coupling graph image showing the connections between the
        OpenFOAM application and FMUs using the Mermaid library.
        """
        from python_mermaid.diagram import MermaidDiagram, Node, Link

        nodes = {
            "openfoam": Node("OpenFOAM"),
            "fmu": Node(self.pyClassName)
        }

        fromPorts = self.get_ports_names("from")
        toPorts = self.get_ports_names("to")

        links  = [Link(nodes["openfoam"], nodes["fmu"], message=port) for port in fromPorts]
        links += [Link(nodes["fmu"], nodes["openfoam"], message=port) for port in toPorts]

        chart = MermaidDiagram(
            title="FMI Coupling",
            nodes=list(nodes.values()),
            links=links
        )

        generate_mermaid_graph_as_image(f"{chart}", "fig_graph_coupling_fmi.png")

        return(chart)



class FunctionObjects(OpenFOAMListDict):
    """
    Object that stores all the functionObjects

    Parameters
    ----------
    functionObjects: list[FunctionObject]
        List of function objects
    """
    def __init__(
            self,
            functionObjects: list[FunctionObject]=None
        ):
        super().__init__(FunctionObject, "functions", functionObjects)
