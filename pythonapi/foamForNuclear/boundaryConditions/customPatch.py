from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.common import Vector


class CustomPatch(Patch):
    """
    CustomPatch patch. Add all the parameters and type of the patch using the
    `parameters` argument.
    """

    def __init__(
            self,
            parameters: dict={},
            value: float | int | Vector = None
        ):
        super().__init__(type="customPatch", value=value)

        for key, item in parameters.items():
            self.__setitem__(key, item)
