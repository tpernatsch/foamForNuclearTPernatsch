from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class NoSlip(Patch):
    """
    This boundary condition provides a no slip constraint.
    """

    def __init__(self):
        super().__init__(type="noSlip")
