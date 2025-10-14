from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict


_DISPERSED_DIAMETER_MODEL_TYPES = {"constant"}


class DispersedDiameterModel(OpenFOAMDict):
    def __init__(self, type: str):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _DISPERSED_DIAMETER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__('type', type)


class ConstantDispersedDiameterModel(DispersedDiameterModel):
    def __init__(
            self,
            value: float
        ):
        super().__init__(type="constant")

        self.value = value

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, value) -> None:
        check_type("value", value, (float, int))
        self._value = value
        self.__setitem__('value', value)
