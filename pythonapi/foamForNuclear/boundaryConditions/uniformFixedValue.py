import numpy as np
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Table
from .boundaryCondition import Patch


class UniformFixedValue(Patch):
    """
    This boundary condition provides a uniform fixed value condition.

    Parameters
    ----------
    uniformValue : Table
    """

    def __init__(self, uniformValue: Table):
        super().__init__(type="uniformFixedValue")

        self.uniformValue: Table = uniformValue


    def __repr__(self, depth=0):
        self.__setitem__("uniformValue", self.uniformValue)

        return super().__repr__(depth)


    @property
    def uniformValue(self):
        return self._uniformValue

    @uniformValue.setter
    def uniformValue(self, uniformValue) -> None:
        check_type("uniformValue", uniformValue, (Table, list, np.ndarray))
        if (isinstance(uniformValue, Table)):
            self._uniformValue = uniformValue
        else:
            self._uniformValue = Table(uniformValue)
