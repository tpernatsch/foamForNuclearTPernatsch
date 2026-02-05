from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class FixedValue(Patch):
    """
    Patch boundary condition that applies a **uniform fixed value**.

    The specified value is imposed uniformly over the entire patch and remains
    constant in time.

    Options
    -------
    value : scalar | Vector
        Value prescribed uniformly on the patch.
        (required: True)
    """
    def __init__(self, value: float):
        super().__init__(type="fixedValue", value=value)
