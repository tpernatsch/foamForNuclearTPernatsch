from foamForNuclear.checkvalue import check_type, check_value
from .boundaryCondition import Patch


class GapPressure(Patch):
    """
    Boundary confition that applies a pressure equal to the gap pressure 
    (i.e. the pressure calculated by the `gapGas` model).

    Options
    -------
    value : list | Vector
        Initial displacement value on the patch.
        (required: True)

    relax : scalar
        Relaxation factor for gradient updates.
        (default: 1; required: False)
    """

    def __init__(
            self,
            value,
            relax: float=1
        ):
        super().__init__(type="gapPressure", value=value)
        
        self.relax = relax

    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (int, float))
        self._relax = relax
        self.__setitem__('relax', relax)
