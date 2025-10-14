from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Vector


class CoupledUniformExternalValue(Patch):
    """
    This boundary condition provides a uniform fixed value condition.

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    inputName : str
        FMI port name
    initValue : float | Vector
        Initial value
    """

    def __init__(self, value, inputName: str, initValue: float=None):
        super().__init__(type="coupledUniformExternalValue", value=value)

        self.inputName = inputName
        self.initValue = initValue


    @property
    def inputName(self):
        return self._inputName

    @inputName.setter
    def inputName(self, inputName) -> None:
        check_type("inputName", inputName, str)
        self._inputName = inputName
        self.__setitem__('inputName', self.inputName)

    @property
    def initValue(self):
        return self._initValue

    @initValue.setter
    def initValue(self, initValue) -> None:
        check_type("initValue", initValue, (int, float, Vector))
        self._initValue = initValue
        self.__setitem__("initValue", initValue)
