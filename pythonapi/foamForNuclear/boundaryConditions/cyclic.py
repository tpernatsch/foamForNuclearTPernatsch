from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class Cyclic(Patch):
    """
    Cyclic patch.
    """

    def __init__(self):
        super().__init__(type="cyclic")
