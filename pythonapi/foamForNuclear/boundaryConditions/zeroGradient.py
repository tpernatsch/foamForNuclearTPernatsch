from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class ZeroGradient(Patch):
    """
    This boundary condition applies a zero-gradient condition from the patch.
    """

    def __init__(self):
        super().__init__(type="zeroGradient")
