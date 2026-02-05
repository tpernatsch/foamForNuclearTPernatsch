from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class Symmetry(Patch):
    """
    Symmetry patch.
    """

    def __init__(self):
        super().__init__(type="symmetry")
