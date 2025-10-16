from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class Calculated(Patch):
    """
    This boundary condition is not designed to be evaluated. It is assmued that
    the value is assigned via field assignment, and not via a call to e.g.
    updateCoeffs or evaluate.

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    """

    def __init__(self, value: float):
        super().__init__(type="calculated", value=value)
