from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class Mixed(Patch):
    """

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    """

    def __init__(
            self,
            value: float,
            refValue: float,
            refGradient: float,
            valueFraction: float,
        ):
        super().__init__(type="mixed", value=value)

        self.refGradient = refGradient
        self.valueFraction = valueFraction
        self.refValue = refValue

    @property
    def refValue(self):
        return self._refValue

    @refValue.setter
    def refValue(self, refValue) -> None:
        check_type("refValue", refValue, (float, int))
        self._refValue = refValue
        self.__setitem__('refValue', f"uniform {refValue}")

    @property
    def refGradient(self):
        return self._refGradient

    @refGradient.setter
    def refGradient(self, refGradient) -> None:
        check_type("refGradient", refGradient, (float, int))
        self._refGradient = refGradient
        self.__setitem__('refGradient', f"uniform {refGradient}")

    @property
    def valueFraction(self):
        return self._valueFraction

    @valueFraction.setter
    def valueFraction(self, valueFraction) -> None:
        check_type("valueFraction", valueFraction, (float, int))
        self._valueFraction = valueFraction
        self.__setitem__('valueFraction', f"uniform {valueFraction}")
