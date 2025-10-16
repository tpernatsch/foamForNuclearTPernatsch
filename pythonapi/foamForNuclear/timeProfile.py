import numpy as np
from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict, Table


_TIME_PROFILE_TYPES = {"table", "fmi"}

_OFFBEAT_TIME_PROFILE_TYPES = {"table"}
_OFFBEAT_OUT_OF_BOUNDS_TYPES = {"clamp"}
_OFFBEAT_INTERPOLATION_SCHEME_TYPES = {"linear"}


class TimeProfile(OpenFOAMDict):
    """
    Time profile object. The first value is maintained constant from -inf to the
    first time and the last value is maintained constant from the last time to
    +inf.

    Example ::

        TableProfile = TimeProfile(
            'table',
            startTime=100,
            table=[(0, 0), (1.1, 57e-5)]
        )
        FmiProfile = TimeProfile(
            'fmi',
            nameFromFMU="gfExtReact",
            initialValue=0
        )

    Parameters
    ----------
    type : str
        Type of time profile (options: `table`, `fmi`)
    startTime : float
        Start time if use type table (default to 0).
    table : list[tuple] | Table
        Table containing the time profile data if use type table
    nameFromFMU : str
        Name of the FMI port
    initialValue : float
        Initial value of the FMI port
    """

    def __init__(
            self,
            type: str,
            table: list | np.ndarray | Table=None,
            startTime: float=0,
            nameFromFMU: str=None,
            initialValue: float=0
        ):
        super().__init__()
        self.type = type
        self.startTime = startTime
        self.table = table
        self.nameFromFMU = nameFromFMU
        self.initialValue = initialValue


    def __repr__(self, depth = 0):
        if (self.type == 'table' and self.table is not None):
            self.__setitem__("startTime", self.startTime)
            self.__setitem__("table", self.table)

        elif (self.type == 'fmi' and self.nameFromFMU is not None):
            self.__setitem__("nameFromFMU", self.nameFromFMU)
            self.__setitem__("initialValue", self.initialValue)

        return super().__repr__(depth)


    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type):
        check_type("type", type, str)
        check_value("type", type, _TIME_PROFILE_TYPES)
        self._type_ = type
        self.__setitem__("type", self.type)

    @property
    def startTime(self):
        return self._startTime

    @startTime.setter
    def startTime(self, startTime):
        check_type("startTime", startTime, (float, int))
        self._startTime = startTime

    @property
    def table(self):
        return self._table

    @table.setter
    def table(self, table):
        check_type("table", table, (list, np.ndarray, Table), none_ok=True)
        if (isinstance(table, Table) or table is None):
            self._table = table
        else:
            self._table = Table(table)

    @property
    def nameFromFMU(self):
        return self._nameFromFMU

    @nameFromFMU.setter
    def nameFromFMU(self, nameFromFMU):
        check_type("nameFromFMU", nameFromFMU, str, none_ok=True)
        self._nameFromFMU = nameFromFMU

    @property
    def initialValue(self):
        return self._initialValue

    @initialValue.setter
    def initialValue(self, initialValue):
        check_type("initialValue", initialValue, (float, int))
        self._initialValue = initialValue

    @property
    def is_empty(self) -> bool:
        if (self.type == "table"):
            return(self.table.is_empty)
        return(False)


class OffbeatTimeProfile(OpenFOAMDict):
    """
    Time profile object

    Example ::

        TableProfile = OffbeatTimeProfile(
            'table',
            values=[(0, 0), (1.1, 57e-5)]
        )

    Parameters
    ----------
    type : str
        Type of time profile (options: `table`)
    values : list[tuple] | Table
        Table containing the time profile data if use type table
    outOfBounds : str
        Out of bound type (options: `clamp`)
    interpolationScheme : str
        Interpolation scheme type (options: `linear`)
    """
    def __init__(
            self,
            type: str,
            values: Table | list[tuple],
            outOfBounds: str="clamp",
            interpolationScheme: str="linear",
        ):
        super().__init__()
        self.type = type
        self.values = values
        self.outOfBounds = outOfBounds
        self.interpolationScheme = interpolationScheme

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type):
        check_type("type", type, str)
        check_value("type", type, _OFFBEAT_TIME_PROFILE_TYPES)
        self._type_ = type
        self.__setitem__("type", self.type)

    @property
    def values(self):
        return self._values

    @values.setter
    def values(self, values):
        check_type("values", values, (Table, list))
        if (isinstance(values, Table)):
            self._values = values
        else:
            self._values = Table(values, type="")
        self.__setitem__("values", self.values)

    @property
    def outOfBounds(self):
        return self._outOfBounds

    @outOfBounds.setter
    def outOfBounds(self, outOfBounds):
        check_type("outOfBounds", outOfBounds, str)
        check_value("outOfBounds", outOfBounds, _OFFBEAT_OUT_OF_BOUNDS_TYPES)
        self._outOfBounds = outOfBounds
        self.__setitem__("outOfBounds", self.outOfBounds)

    @property
    def interpolationScheme(self):
        return self._interpolationScheme

    @interpolationScheme.setter
    def interpolationScheme(self, interpolationScheme):
        check_type("interpolationScheme", interpolationScheme, str)
        check_value("interpolationScheme", interpolationScheme, _OFFBEAT_INTERPOLATION_SCHEME_TYPES)
        self._interpolationScheme = interpolationScheme
        self.__setitem__("interpolationScheme", self.interpolationScheme)
