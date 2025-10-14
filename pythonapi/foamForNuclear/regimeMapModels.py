from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import List, OpenFOAMDict


_REGIME_MAP_MODEL_TYPES = {"oneParameter"}
_INTERPOLATION_MODE_TYPES = {"linear", "quadratic"}


class RegimeMapModel(OpenFOAMDict):
    """
    Base class for regime map model.
    """
    def __init__(self, name: str, type: str):
        super().__init__(name=name)
        self.type = type


    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _REGIME_MAP_MODEL_TYPES)
        self._type_ = type
        self.__setitem__("type", type)


class RegimeMapOneParameter(RegimeMapModel):
    """
    Regime map model dependent on one parameter.
    """
    def __init__(
            self,
            name: str,
            parameter: str,
            regimeBounds: OpenFOAMDict=None,
            interpolationMode: str=None
        ):
        super().__init__(name=name, type="oneParameter")

        self.parameter = parameter
        self.regimeBounds = regimeBounds
        self.interpolationMode = interpolationMode


    def __repr__(self, depth = 0):
        self.__setitem__("regimeBounds", self.regimeBounds)

        return super().__repr__(depth)


    @property
    def parameter(self):
        return self._parameter

    @parameter.setter
    def parameter(self, parameter) -> None:
        check_type("parameter", parameter, str)
        self._parameter = parameter
        self.__setitem__("parameter", parameter)


    @property
    def interpolationMode(self):
        return self._interpolationMode

    @interpolationMode.setter
    def interpolationMode(self, interpolationMode) -> None:
        check_type("interpolationMode", interpolationMode, str, none_ok=True)
        self._interpolationMode = interpolationMode
        if (interpolationMode is not None):
            check_value("interpolationMode", interpolationMode, _INTERPOLATION_MODE_TYPES)
            self.__setitem__("interpolationMode", interpolationMode)


    @property
    def regimeBounds(self):
        return self._regimeBounds

    @regimeBounds.setter
    def regimeBounds(self, regimeBounds) -> None:
        check_type("regimeBounds", regimeBounds, OpenFOAMDict, none_ok=True)
        if (regimeBounds is not None):
            self._regimeBounds = regimeBounds
        else:
            self._regimeBounds = OpenFOAMDict()


    def add_regime(self, name: str, minValue: float, maxValue: float):
        check_type("name", name, str)
        check_type("minValue", minValue, (float, int))
        check_type("maxValue", maxValue, (float, int))
        self.regimeBounds[name] = List([minValue, maxValue])
