from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class FixedFluxExtrapolatedPressure(Patch):
    """
    This boundary condition sets the pressure gradient to the provided value
    such that the flux on the boundary is that specified by the velocity
    boundary condition.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    """

    def __init__(self, value: float, gradient: float):
        super().__init__(type="fixedFluxExtrapolatedPressure", value=value)

        self.gradient = gradient

    @property
    def gradient(self):
        return self._gradient

    @gradient.setter
    def gradient(self, gradient) -> None:
        check_type("gradient", gradient, (float, int))
        self._gradient = gradient
        self.__setitem__('gradient', f"uniform {gradient}")
