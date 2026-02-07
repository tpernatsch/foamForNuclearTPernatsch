import pandas as pd
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import *
import attrs as attr
from attrs import define
from foamForNuclear._attrs_tools import (
    auto_type_validator, _to_List_str, mirror_to_dict, _to_Vector)
import os
import numpy as np
import re

# Duplicated from controlDict.py
_WRITE_CONTROL_TYPES = {
    'none', 'timeStep', 'runTime', 'adjustable', 'adjustableRunTime',
    'writeTime', 'clockTime', 'cpuTime'
}


_FUNCTION_OBJECT_TYPES = {
    'fieldMinMax', 'surfaceFieldValue', 'volFieldValue',
    'mag',
    'probes',
    'patchProbes',
    'sets',
    'multiply', 'subtract', 'divide', 'limitFields',
    'massFlow', 'TBulk',
    'FMUSimulator',
    'fgr'
}
_FUNCTION_OBJECT_LIBS = {
    "fieldFunctionObjects", "libfieldFunctionObjects.so", "libsampling.so",
    'pyFMUSim', "libFunctionObjects.so", "libFFNFunctionObjects.so", "libOffbeatFunctionObject.so"
}
_VOLUME_OPERATION_TYPES = {
    "none", "min", "max", "sum", "sumMag", "average", "volAverage", 
    "volIntegrate", "CoV", "weightedSum", "weightedAverage", 
    "weightedVolAverage", "weightedVolIntegrate"
}
_VOLUME_REGION_TYPES = {
    "cellZone", "all"
}
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

_MODE_TYPES = {"magnitude", "component"}

_PAREN_GROUP = re.compile(r"\(([^()]*)\)")

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
            scaleFactor: float=None,
            caseFolder: str=None
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
        self.caseFolder = caseFolder


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

    @property
    def caseFolder(self):
        return self._caseFolder

    @caseFolder.setter
    def caseFolder(self, caseFolder):
        check_type("caseFolder", caseFolder, str, none_ok=True)
        self._caseFolder = caseFolder
        if (caseFolder is not None):
            self.__setitem__("caseFolder", caseFolder)

    def read_from_case(self, startTime: float, caseFolder: str | None =None):
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
        OpenFOAM's Maps input fields from local mesh to secondary mesh at
        runtime.

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


class FieldMinMax(FunctionObject):
    def __init__(
            self,
            name: str,
            fields: list[str],
            mode: str,
            location: bool=True,
            log: bool=None,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
        ):
        """
        The OpenFOAM's fieldMinMax function object computes the values and
        locations of field minima and maxima. These are good indicators of
        calculation performance, e.g. to confirm that predicted results are
        within expected bounds, or how well a case is converging.

        Multiple fields can be processed, where for rank > 0 primitives, e.g.
        vectors and tensors, the extrema can be calculated per component, or by
        magnitude. In addition, spatial location and local processor index are
        included in the output.

        More info: https://www.openfoam.com/documentation/guides/latest/doc/guide-fos-field-fieldMinMax.html

        Parameters
        ----------
        fields : list[str]
            List of fields name to map from `region` to `mapRegion`.
        mode : str {magnitude, component}
            Calculation mode: magnitude or component.
        location : bool
            Write location of the min/max value.
        region : str
            Name of the region to map from.
        """
        super().__init__(
            name,
            "fieldMinMax",
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
        self.mode = mode
        self.location = location

    @property
    def fields(self):
        return self._fields

    @fields.setter
    def fields(self, fields):
        check_type("fields", fields, list)
        self._fields = fields
        self.__setitem__("fields", List(fields))

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, mode):
        check_type("mode", mode, str)
        check_value("mode", mode, _MODE_TYPES)
        self._mode = mode
        self.__setitem__("mode", mode)

    @property
    def location(self):
        return self._location

    @location.setter
    def location(self, location):
        check_type("location", location, bool)
        self._location = location
        self.__setitem__("location", location)

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
            caseFolder: str | None = None
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time

        caseFolder : str | None
            Case folder to read from (default None; read from own caseFolder attribute).
        """
        if caseFolder is None:
            caseFolder = getattr(self, "caseFolder", None)
        if caseFolder is None:
            raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")
        
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
            writeFields: bool=False,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
            regionType: str="cellZone",
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
    def regionName(self):
        return self._regionName

    @regionName.setter
    def regionName(self, regionName):
        check_type("regionName", regionName, str)
        self._regionName = regionName
        self.__setitem__("name", regionName)

    @property
    def operation(self):
        return self._operation

    @operation.setter
    def operation(self, operation):
        check_type("operation", operation, str)
        check_value("operation", operation, _VOLUME_OPERATION_TYPES)
        self._operation = operation
        self.__setitem__("operation", operation)    

    @property
    def regionType(self):
        return self._regionType

    @regionType.setter
    def regionType(self, regionType):
        check_type("regionType", regionType, str)
        check_value("regionType", regionType, _VOLUME_REGION_TYPES)
        self._regionType = regionType
        self.__setitem__("regionType", regionType)

    def read_from_case(
            self,
            startTime: float,
            caseFolder: str | None = None
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time (folder name under postProcessing).
        caseFolder : str | None
            Case folder to read from (default None; read from own caseFolder attribute).

        Returns
        -------
        times : list[float]
            List of time values.
        data : dict[str, list]
            Keys are clean field names (e.g. 'D', 'Bu', 'T'):
            - scalar  -> list[float]
            - vector/tensor -> list[list[float]]
        """
        if caseFolder is None:
            caseFolder = getattr(self, "caseFolder", None)
        if caseFolder is None:
            raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")
        
        region = "" if self.region is None else self.region

        folder = os.path.join(
            caseFolder,
            "postProcessing",
            region,
            self.name,
            str(startTime)
        )
        filename = os.path.join(folder, "volFieldValue.dat")

        times: list[float] = []
        data: dict[str, list] = {}

        def _clean_name(raw: str) -> str:
            # e.g. "volAverage(Bu)" -> "Bu"
            if raw.endswith(")") and "(" in raw:
                return raw[raw.find("(") + 1 : -1]
            return raw

        with open(filename, "r") as f:
            header_cols = None

            # --- find header line ("# Time ...") ---
            for line in f:
                if not line.strip():
                    continue
                if line.startswith("# Time"):
                    header_line = line.lstrip("#").strip()
                    header_cols = header_line.split()
                    # header_cols[0] = "Time", others are e.g. "volAverage(D)"
                    raw_fields = header_cols[1:]
                    field_names = [_clean_name(name) for name in raw_fields]
                    data = {name: [] for name in field_names}
                    break

            if header_cols is None:
                raise RuntimeError(f"No header line '# Time ...' found in {filename}")

            # --- read data lines ---
            for line in f:
                if not line.strip():
                    continue

                parts = line.strip().split("\t")
                if len(parts) < 2:
                    continue

                # time column
                t = float(parts[0].strip())
                times.append(t)

                # other columns: one per field
                for name, val_str in zip(field_names, parts[1:]):
                    v = val_str.strip()

                    # vector / tensor: "(x y z ...)"
                    if v.startswith("(") and v.endswith(")"):
                        inner = v[1:-1].strip()
                        comps = [float(c) for c in inner.split()]
                        data[name].append(comps)
                    else:
                        # scalar
                        data[name].append(float(v))

        return data, times



class Probes(FunctionObject):
    """
    Probes functionObject
    """
    def __init__(
            self,
            name,
            fields: list[str],
            enabled: bool,
            probeLocations: list[list[float]] | list[tuple] | list[Vector],
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

        # User-facing container: behaves like a list, type-checked
        self.probeLocations: list[list[float]] | list[tuple] | list[Vector] = CheckedList(
            (list, tuple, Vector, np.ndarray), "probeLocations"
        )

        # Initialise from constructor argument via the helper (for back-compat)
        if probeLocations is not None:
            for pl in probeLocations:
                self.add_probe_location(pl)

    def __repr__(self, depth=0):
        # Build a *separate* list of Vectors for OpenFOAM representation
        vec_locs = []
        for pl in self.probeLocations:
            if isinstance(pl, Vector):
                vec_locs.append(pl)
            else:
                # tuple / list / ndarray → Vector
                v = Vector(pl[0], pl[1], pl[2])
                vec_locs.append(v)

        # Use your List wrapper just for the dict value
        self.__setitem__("probeLocations", List(vec_locs))
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
            caseFolder: str | None = None
        ):
        if caseFolder is None:
            caseFolder = getattr(self, "caseFolder", None)
        if caseFolder is None:
            raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")

        region = "" if self.region is None else self.region
        filename = f"{caseFolder}/postProcessing/{self.name}/{region}/{startTime}/{fieldName}"

        locations = []
        data = {}

        with open(filename, "r") as f:
            for line in f:
                s = line.strip()
                if not s:
                    continue

                if "# Probe " in s:
                    pos = s.replace("(", "").replace(")", "").split()[-3:]
                    x, y, z = (float(v) for v in pos)
                    locations.append(Vector(x, y, z))
                    continue

                if s.startswith("#"):
                    continue

                # ---- data line ----
                first, *rest = s.split(maxsplit=1)
                time = float(first)
                rest = rest[0] if rest else ""

                groups = _PAREN_GROUP.findall(rest)

                if groups:
                    # vector / tensor
                    vals = np.array(
                        [[float(v) for v in g.split()] for g in groups],
                        dtype=float
                    )
                else:
                    # scalar
                    vals = np.array(
                        [float(v) for v in rest.split()] if rest else [],
                        dtype=float
                    )

                # squeeze probe dimension if single probe
                if vals.ndim > 1 and vals.shape[0] == 1:
                    vals = vals[0]

                data[time] = vals

        return data, locations

   
class PatchProbes(Probes):
    """
    PatchProbes functionObject
    """
    def __init__(
            self,
            name,
            fields: list[str],
            enabled: bool,
            probeLocations: list[list[float]] | list[tuple] | list[Vector],
            patches: list[str],
            log: bool=None,
            writeFields: bool=None,
            writeControl: str=None,
            writeInterval: float=None,
            region: str=None,
        ):
        super().__init__(
            name,
            fields,
            enabled,
            probeLocations,
            log,
            writeFields,
            writeControl,
            writeInterval,
            region,
        )
        self.type = "patchProbes"

        # User-facing container: behaves like a list, type-checked
        self.patches: list[str] = CheckedList(
            (str), "patches"
        )

        # Initialise from constructor argument via the helper (for back-compat)
        if patches is not None:
            for p in patches:
                self.patches.append(p)

    def __repr__(self, depth=0):
        of_version = of_flavour()
        # Use your List wrapper just for the dict value
        if(of_version == "esi"):
            self.__setitem__("patches", List(self.patches))

        elif(of_version == "foundation"):
            if len(self.patches) != 1:
                raise ValueError(
                    "OpenFOAM Foundation patchProbes supports only one patch name. "
                    f"Got {len(self.patches)} patches: {self.patches}"
                )
            self.__setitem__("patchName", str(self.patches[0]))

        return super().__repr__(depth)

    @property
    def patches(self):
        return self._patches

    @patches.setter
    def patches(self, patches):
        check_type("patches", patches, list)
        self._patches = patches
        self.__setitem__("patches", List(patches))

    # def read_from_case(
    #         self,
    #         startTime: float,
    #         fieldName: str,
    #         caseFolder: str | None =None
    #     ):
    #     """
    #     Parameters
    #     ----------
    #     startTime : float
    #         Simulation start time
    #     fieldName : str
    #         Name of the field
    #     caseFolder : str | None
    #         Case folder to read from (default None; read from own caseFolder attribute).

    #     Return
    #     ------
    #     Return the probe values and the data locations
    #     """
    #     if caseFolder is None:
    #         caseFolder = getattr(self, "caseFolder", None)
    #     if caseFolder is None:
    #         raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")
        
    #     region = "" if self.region is None else self.region
    #     filename = f'{caseFolder}/postProcessing/{self.name}/{region}/{startTime}/{fieldName}'

    #     locations = []
    #     data = {}
    #     with open(filename, 'r') as f:
    #         for line in f:
    #             if ("# Probe " in line):
    #                 pos = line.replace('(', '').replace(')', '').split()[-3:]
    #                 pos = [float(val) for val in pos]
    #                 vec = Vector(pos[0], pos[1], pos[2])
    #                 locations.append(vec)
    #             elif ("# Time " in line):
    #                 continue
    #             else:
    #                 row = [float(val) for val in line.split()]
    #                 time, values = row[0], row[1:]
    #                 data[time] = np.array(values)

    #     return(data, locations)


@define(slots=True, on_setattr=[attr.setters.convert, attr.setters.validate],
        field_transformer=auto_type_validator, repr=False, kw_only=True)
class Graph(FunctionObject):
    """
    Graph functionObject
    """
    name: str
    start: tuple | list | Vector = attr.field(converter=_to_Vector)
    end: tuple | list | Vector = attr.field(converter=_to_Vector)
    fields: list = attr.field(converter=_to_List_str)
    graph_type: str = "uniform"
    axis: str = "distance"
    nPoints: float | int = 100
    interpolationScheme: str = "cellPoint"
    writeControl: str = "writeTime"
    writeInterval: str | None = None
    setFormat: str = "csv"
    writeFields: bool | None = None
    log: bool | None = None
    region: str | None = None

    def __attrs_post_init__(self):
        super().__init__(
            self.name,
            "sets",
            "libsampling.so",
            self.log,
            self.writeFields,
            self.writeControl,
            self.writeInterval,
            self.region,
            regionType=None,
            regionName=None,
            scaleFactor=None
        )

    def __repr__(self, depth=0):
        of_version = of_flavour()
        
        if(of_version == "foundation" and self.graph_type == "uniform"):
            self.graph_type = "lineUniform"

        setConfig = OpenFOAMDict({})
        line = OpenFOAMDict(
            name='line', items={'type': self.graph_type, 'axis': self.axis, 'nPoints': self.nPoints, 'start': self.start, 'end': self.end})
        sets = OpenFOAMList(name='sets', expected_type=OpenFOAMDict,items=[line])

        self.__setitem__('fields', self.fields)
        self.__setitem__('writeControl', self.writeControl)
        self.__setitem__('interpolationScheme',  self.interpolationScheme)
        self.__setitem__('setFormat',  self.setFormat)
        self.__setitem__('sets', sets)

        return super().__repr__(depth)

    def read_from_case(self, startTime: float, caseFolder: str | None = None):
        if caseFolder is None:
            caseFolder = getattr(self, "caseFolder", None)
        if caseFolder is None:
            raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")

        base = os.path.join(caseFolder, "postProcessing", self.name)
        if not os.path.isdir(base):
            raise RuntimeError(f"Post-processing folder not found: {base}")

        data = {fld: {} for fld in self.fields}
        locations = None

        def _read_one_file(fpath: str):
            with open(fpath, "r") as f:
                header_line = f.readline().strip()
            header = [h.strip().strip('"') for h in header_line.split(",")]

            arr = np.loadtxt(fpath, delimiter=",", skiprows=1)
            if arr.ndim == 1:
                arr = arr.reshape(1, -1)

            x = arr[:, 0]
            return header, arr, x

        def _extract_field_from_header_arr(header, arr, field: str):
            cols = []
            for i, h in enumerate(header):
                if i == 0:
                    continue
                if h == field or h.startswith(field + "_"):
                    cols.append(i)
            if not cols:
                return None

            vals = arr[:, cols]
            if vals.ndim == 1 or vals.shape[1] == 1:
                return vals.flatten().tolist()
            return [row.tolist() for row in vals]

        # collect times >= startTime
        times = []
        for d in os.listdir(base):
            try:
                t = float(d)
            except ValueError:
                continue
            if t >= startTime:
                times.append((t, d))
        times.sort()

        for t, tdir in times:
            time_path = os.path.join(base, tdir)

            # ---------------------------
            # 1) Try ESI-style: single file containing all fields (encoded in filename)
            # ---------------------------
            candidates = []
            for fname in os.listdir(time_path):
                if not (fname.startswith("line_") and fname.endswith(".csv")):
                    continue

                base_name = fname[:-4]      # strip .csv
                name_body = base_name[5:]   # strip "line_"

                ok = True
                for fld in self.fields:
                    if not (
                        name_body == fld
                        or name_body.startswith(fld + "_")
                        or name_body.endswith("_" + fld)
                        or ("_" + fld + "_") in name_body
                    ):
                        ok = False
                        break
                if ok:
                    candidates.append(fname)

            used_mode = None
            if len(candidates) == 1:
                used_mode = "single"
                fpath = os.path.join(time_path, candidates[0])

                header, arr, x = _read_one_file(fpath)

                if locations is None:
                    locations = x.tolist()
                else:
                    if not np.allclose(locations, x):
                        raise RuntimeError("Inconsistent radial coordinates across times")

                for field in self.fields:
                    values_list = _extract_field_from_header_arr(header, arr, field)
                    if values_list is not None:
                        data[field][t] = values_list

            elif len(candidates) > 1:
                raise RuntimeError(
                    f"More than one line file matches fields {self.fields} "
                    f"at time {t}: {candidates}"
                )

            # ---------------------------
            # 2) If no single file, try Foundation-style: one file per field: line_<field>.csv
            # ---------------------------
            if used_mode is None:
                any_found = False
                for field in self.fields:
                    fname = f"line_{field}.csv"
                    fpath = os.path.join(time_path, fname)
                    if not os.path.isfile(fpath):
                        continue

                    any_found = True
                    header, arr, x = _read_one_file(fpath)

                    if locations is None:
                        locations = x.tolist()
                    else:
                        if not np.allclose(locations, x):
                            raise RuntimeError("Inconsistent radial coordinates across times")

                    values_list = _extract_field_from_header_arr(header, arr, field)
                    if values_list is None:
                        raise RuntimeError(
                            f"Field '{field}' not found in header of {fname}: {header}"
                        )
                    data[field][t] = values_list

                # If neither mode found anything for this time, just skip (consistent with old behavior)
                if not any_found:
                    continue

        if locations is None:
            raise RuntimeError(
                f"No matching line files found in {base} for times >= {startTime}"
            )

        return data, locations




@define(
    slots=True,
    on_setattr=[attr.setters.convert, attr.setters.validate],
    field_transformer=auto_type_validator,
    repr=False,
    kw_only=True,
)
class FGR(FunctionObject):
    """
    fgr functionObject

    Example output:
        fgr
        {
            type fgr;
        }
    """
    name: str

    # keep the usual FO “meta” options optional
    writeControl: str | None = None
    writeInterval: str | None = None
    writeFields: bool | None = None
    log: bool | None = None
    region: str | None = None

    def __attrs_post_init__(self):
        super().__init__(
            self.name,
            "fgr",
            "libOffbeatFunctionObject.so",               # no library specified (matches your minimal dict)
            self.log,
            self.writeFields,
            self.writeControl,
            self.writeInterval,
            self.region,
            regionType=None,
            regionName=None,
            scaleFactor=None,
        )

    def read_from_case(
            self,
            startTime: float,
            caseFolder: str | None = None,
            filename: str = "fgr",
        ):
        """
        Read FGR time series from:
            postProcessing/<self.name>/<region>/<startTime>/<filename>

        File format:
            # Time fgr(%)
            t0  v0
            t1  v1
            ...

        Returns
        -------
        data : dict[float, float]
            data[time] = fgr_value
        """
        if caseFolder is None:
            caseFolder = getattr(self, "caseFolder", None)
        if caseFolder is None:
            raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")

        region = "" if self.region is None else self.region
        fpath = f"{caseFolder}/postProcessing/{self.name}/{region}/{startTime}/{filename}.dat"

        data = {}

        with open(fpath, "r") as f:
            for line in f:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue

                parts = s.split()
                if len(parts) < 2:
                    continue

                t = float(parts[0])
                v = float(parts[1])
                data[t] = v

        return data



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
            caseFolder: str | None = None
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time
        caseFolder : str | None
            Case folder to read from (default None; read from own caseFolder attribute).
        """
        if caseFolder is None:
            caseFolder = getattr(self, "caseFolder", None)
        if caseFolder is None:
            raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")
        
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
            caseFolder: str | None = None
        ):
        """
        Parameters
        ----------
        startTime : float
            Simulation start time
        caseFolder : str | None
            Case folder to read from (default None; read from own caseFolder attribute).
        """
        if caseFolder is None:
            caseFolder = getattr(self, "caseFolder", None)
        if caseFolder is None:
            raise ValueError("Need caseFolder or a bound caseFolder on this functionObject.")
        
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
    owner: OffbeatCase | None
        Parent case, used to stamp caseFolder on functions
    """
    def __init__(
            self,
            functionObjects: list[FunctionObject] = None,
            owner='./',
        ):
        # Initialise with an *empty* list; we'll extend manually to bind owner caseFolder
        super().__init__(FunctionObject, "functions", [])

        self.owner = owner

        if functionObjects:
            self.extend(functionObjects)

    # ---- internal helper to bind owner caseFolder ----

    def _bind_case_folder(self, func: FunctionObject) -> None:
        """Stamp the owner's caseFolder on the function, if available."""
        if self.owner is None:
            return
        folder = getattr(self.owner, "caseFolder", None)
        if folder is None:
            return
        # Just set an attribute; FunctionObject.read_from_case will use it
        setattr(func, "caseFolder", folder)

    # ---- overrides that add/replace items ----

    def append(self, item: FunctionObject) -> None:
        super().append(item)
        # bind the last stored element, not the raw argument
        self._bind_case_folder(self[-1])

    def extend(self, iterable) -> None:
        start = len(self)
        super().extend(iterable)
        # bind all *newly added* items from the container itself
        for f in self[start:]:
            self._bind_case_folder(f)


    def insert(self, index: int, item: FunctionObject) -> None:
        super().insert(index, item)
        self._bind_case_folder(self[index])

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if isinstance(key, slice):
            # bind everything now stored in that slice
            for f in self[key]:
                self._bind_case_folder(f)
        else:
            self._bind_case_folder(self[key])


    def __iadd__(self, other):
        start = len(self)
        super().__iadd__(other)
        for f in self[start:]:
            self._bind_case_folder(f)
        return self


    # ---- export ----

    def export_to_openfoam(self, settings):
        if not self:
            return

        # Read current controlDict
        with open(settings.path, "r+", encoding="utf-8") as f:
            content = f.read()

            footer = openfoamFooterLine.strip()
            stripped = content.rstrip()
            had_footer = stripped.splitlines()[-1].strip() == footer if stripped else False

            if had_footer:
                # Remove the footer (and trailing whitespace) before appending
                head = stripped.rsplit(footer, 1)[0].rstrip()
            else:
                head = stripped

            # Build new content: head + functions + optional footer
            fn_block = repr(self).rstrip()
            new_text = head + "\n\n" + fn_block + "\n"
            if had_footer:
                new_text += "\n" + openfoamFooterLine

            # Overwrite file
            f.seek(0)
            f.write(new_text)
            f.truncate()