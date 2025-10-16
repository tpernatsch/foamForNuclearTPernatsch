from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type


class FixedJump(Patch):
    """
    The fixedJump is a general coupled boundary condition to provide a jump
    cyclic condition. The jump is specified as a fixed value field, applied as
    an offset to the 'owner' patch.

    Parameters
    ----------
    value : float | int
        Value at the patch set uniformly
    jump : float | int
        Jump field
    jump0 : float | int
        Old-time level jump field (default `None`)
    relax : float | int
        Under-relaxation factor (default `1`)
    rho : str
        Density field name (default `thermo:rho`)
    """
    def __init__(
            self,
            value: float,
            jump: float,
            jump0: float=None,
            relax: float=None,
            rho: str="thermo:rho",
            patchType: str="cyclic"
        ):
        super().__init__(type="fixedJump", value=value)

        self.patchType = patchType
        self.jump = jump
        self.jump0 = jump0
        self.relax = relax
        self.rho = rho

    @property
    def jump(self):
        return self._jump

    @jump.setter
    def jump(self, jump) -> None:
        check_type("jump", jump, (float, int))
        self._jump = jump
        self.__setitem__('jump', f"uniform {jump}")

    @property
    def jump0(self):
        return self._jump0

    @jump0.setter
    def jump0(self, jump0) -> None:
        check_type("jump0", jump0, (float, int), none_ok=True)
        self._jump0 = jump0
        if (jump0 is not None):
            self.__setitem__('jump0', f"uniform {jump0}")

    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (float, int), none_ok=True)
        self._relax = relax
        if (relax is not None):
            self.__setitem__('relax', relax)

    @property
    def patchType(self):
        return self._patchType

    @patchType.setter
    def patchType(self, patchType) -> None:
        check_type("patchType", patchType, str)
        self._patchType = patchType
        self.__setitem__('patchType', patchType)

    @property
    def rho(self):
        return self._rho

    @rho.setter
    def rho(self, rho) -> None:
        check_type("rho", rho, str)
        self._rho = rho
        self.__setitem__('rho', rho)
