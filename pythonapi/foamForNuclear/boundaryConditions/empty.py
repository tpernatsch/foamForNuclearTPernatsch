from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class Empty(Patch):
    """
    Empty patch.
    """

    def __init__(self):
        super().__init__(type="empty")
