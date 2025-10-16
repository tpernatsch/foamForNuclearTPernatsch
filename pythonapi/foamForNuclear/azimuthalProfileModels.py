from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import OpenFOAMDict


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
