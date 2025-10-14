from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class Slip(Patch):
    """
    This boundary condition provides a slip constraint.
    """

    def __init__(self):
        super().__init__(type="slip")
