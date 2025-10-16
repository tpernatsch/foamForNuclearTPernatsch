from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class FixedValue(Patch):
    """
    Patch with a uniform distribution.

    Parameters
    ----------
    value : float | int | Vector
        Value at the patch set uniformly
    """

    def __init__(self, value: float):
        super().__init__(type="fixedValue", value=value)
