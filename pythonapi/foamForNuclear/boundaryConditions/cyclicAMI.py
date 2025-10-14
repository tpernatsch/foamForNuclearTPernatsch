from foamForNuclear.boundaryConditions.boundaryCondition import Patch


class CyclicAMI(Patch):
    """
    Cyclic AMI patch.
    """
    def __init__(self):
        super().__init__(type="cyclicAMI")
