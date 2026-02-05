from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import List, OpenFOAMDict, Table

_INTERPOLATION_METHOD_TYPES = {"linear"}


class AxialProfile(OpenFOAMDict):
    """
    Base class axial profile.
    """
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        self._type = type
        self.__setitem__('type', self.type)


class FlatAxialProfile(AxialProfile):
    """
    Flat axial profile.
    """
    def __init__(self):
        super().__init__(type="flat")


class TimeDependentTabulatedAxialProfile(AxialProfile):
    def __init__(
            self,
            timePoints: list | List,
            axialLocations: list | List,
            data: list | Table,
            axialInterpolationMethod: str="linear",
            burnupInterpolationMethod: str="linear",
        ):
        super().__init__(type="timeDependentTabulated")

        self.timePoints = timePoints
        self.axialLocations = axialLocations
        self.data = data
        self.axialInterpolationMethod = axialInterpolationMethod
        self.burnupInterpolationMethod = burnupInterpolationMethod


    def __repr__(self, depth = 0):
        self.__setitem__('timePoints', self.timePoints)
        self.__setitem__('axialLocations', self.axialLocations)
        self.__setitem__('data', self.data)

        return super().__repr__(depth)


    @property
    def timePoints(self):
        return self._timePoints

    @timePoints.setter
    def timePoints(self, timePoints) -> None:
        check_type("timePoints", timePoints, (list, List))
        if (isinstance(timePoints, list)):
            self._timePoints = List(timePoints)
        else:
            self._timePoints = timePoints

    @property
    def axialLocations(self):
        return self._axialLocations

    @axialLocations.setter
    def axialLocations(self, axialLocations) -> None:
        check_type("axialLocations", axialLocations, (list, List))
        if (isinstance(axialLocations, list)):
            self._axialLocations = List(axialLocations)
        else:
            self._axialLocations = axialLocations

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data) -> None:
        check_type("data", data, (list, Table))
        if (isinstance(data, list)):
            self._data = Table(data, type="")
        else:
            self._data = data
            self._data.type = ""

    @property
    def axialInterpolationMethod(self):
        return self._axialInterpolationMethod

    @axialInterpolationMethod.setter
    def axialInterpolationMethod(self, axialInterpolationMethod) -> None:
        check_type("axialInterpolationMethod", axialInterpolationMethod, str)
        check_value("axialInterpolationMethod", axialInterpolationMethod, _INTERPOLATION_METHOD_TYPES)
        self._axialInterpolationMethod = axialInterpolationMethod
        self.__setitem__('axialInterpolationMethod', self.axialInterpolationMethod)

    @property
    def burnupInterpolationMethod(self):
        return self._burnupInterpolationMethod

    @burnupInterpolationMethod.setter
    def burnupInterpolationMethod(self, burnupInterpolationMethod) -> None:
        check_type("burnupInterpolationMethod", burnupInterpolationMethod, str)
        check_value("burnupInterpolationMethod", burnupInterpolationMethod, _INTERPOLATION_METHOD_TYPES)
        self._burnupInterpolationMethod = burnupInterpolationMethod
        self.__setitem__('burnupInterpolationMethod', self.burnupInterpolationMethod)


class AzimuthalProfile(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        self._type = type
        self.__setitem__('type', self.type)


class RadialProfile(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        self._type = type
        self.__setitem__('type', self.type)


class FromBurnupRadialProfile(RadialProfile):
    def __init__(
            self
        ):
        super().__init__(type="fromBurnup")
