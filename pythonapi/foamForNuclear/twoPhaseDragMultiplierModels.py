from foamForNuclear.checkvalue import check_type, check_value
from foamForNuclear.common import OpenFOAMDict


_TWO_PHASE_DRAG_MULTIPLIER_MODEL_TYPES = {"LottesFlinn", "Kaiser88", "LockhartMartinelli"}


class TwoPhaseDragMultiplierModel(OpenFOAMDict):
    def __init__(
            self,
            type: str
        ):
        super().__init__()

        self.type = type

    @property
    def type(self):
        return self._type_

    @type.setter
    def type(self, type) -> None:
        check_type("type", type, str)
        check_value("type", type, _TWO_PHASE_DRAG_MULTIPLIER_MODEL_TYPES)
        self._type_ = type
        self.__setitem__('type', type)


class LottesFlinn(TwoPhaseDragMultiplierModel):
    def __init__(
            self,
            multiplierFluid: str
        ):
        super().__init__("LottesFlinn")

        self.multiplierFluid = multiplierFluid

    @property
    def multiplierFluid(self):
        return self._multiplierFluid

    @multiplierFluid.setter
    def multiplierFluid(self, multiplierFluid) -> None:
        check_type("multiplierFluid", multiplierFluid, str)
        self._multiplierFluid = multiplierFluid
        self.__setitem__('multiplierFluid', multiplierFluid)


class LockhartMartinelli(TwoPhaseDragMultiplierModel):
    def __init__(
            self,
            multiplierFluid: str
        ):
        super().__init__("LockhartMartinelli")

        self.multiplierFluid = multiplierFluid

    @property
    def multiplierFluid(self):
        return self._multiplierFluid

    @multiplierFluid.setter
    def multiplierFluid(self, multiplierFluid) -> None:
        check_type("multiplierFluid", multiplierFluid, str)
        self._multiplierFluid = multiplierFluid
        self.__setitem__('multiplierFluid', multiplierFluid)


class Kaiser88(TwoPhaseDragMultiplierModel):
    def __init__(
            self,
            multiplierFluid: str
        ):
        super().__init__("Kaiser88")

        self.multiplierFluid = multiplierFluid

    @property
    def multiplierFluid(self):
        return self._multiplierFluid

    @multiplierFluid.setter
    def multiplierFluid(self, multiplierFluid) -> None:
        check_type("multiplierFluid", multiplierFluid, str)
        self._multiplierFluid = multiplierFluid
        self.__setitem__('multiplierFluid', multiplierFluid)
