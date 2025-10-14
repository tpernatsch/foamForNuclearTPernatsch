from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class Wedge(Patch):
    """
    This boundary condition is similar to the cyclic condition, except that it
    is applied to 2-D geometries.
    """

    def __init__(self):
        super().__init__(type="wedge")
