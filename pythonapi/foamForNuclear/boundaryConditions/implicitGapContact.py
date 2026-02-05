from foamForNuclear.checkvalue import check_type, check_value
from .boundaryCondition import Patch

_REGION_COUPLED_TYPES = {"regionCoupledOFFBEAT"}

class ImplicitGapContact(Patch):
    """
    Coupled contact boundary condition for displacement, using an **implicit**
    penalty formulation.

    The `ImplicitGapContact` fvPatchField is selected in the patch subdictionary
    inside the `boundaryField` of the displacement field. It is intended to be
    used on a *pair* of `regionCoupledOFFBEAT` patches coupled via AMI.

    The normal contact pressure is enforced with a penalty method. Optional friction
    can be activated via a Coulomb friction coefficient and a separate penalty scaling
    for tangential forces. A constant offset can be applied to the computed gap width.

    Options
    -------
    value : list | Vector
        Initial displacement value on the patch.
        (required: True)

    patchType : word
        Type of the underlying fvPatch. Must be set to `regionCoupledOFFBEAT`.
        (default: regionCoupledOFFBEAT; required: False)

    penaltyFactor : scalar
        Penalty scale factor for the normal contact pressure. Controls the balance
        between penetration and convergence. Larger values correspond to stiffer
        contacts.
        Typical values: 0.01 to 1 (or larger for stiff contacts).
        (default: 0.1; required: False)

    penaltyFactorFriction : scalar
        Penalty scale factor for the frictional (tangential) forces.
        Typical values: 0.01 to 1.
        (default: 0.1; required: False)

    frictionCoefficient : scalar
        Coulomb friction coefficient between surfaces.
        (default: 0; required: False)

    offset : scalar
        Offset added to the computed gap width. A positive value increases the gap.
        (default: 0; required: False)

    glued : bool
        If true, treat the two surfaces as glued (attached), preventing relative
        motion across the interface.
        (default: False; required: False)

    relax : scalar
        Relaxation factor for gradient updates.
        (default: 1; required: False)
    """

    def __init__(
            self,
            value,
            patchType: str="regionCoupledOFFBEAT",
            penaltyFactor: float=0.1,
            penaltyFactorFriction: float=0.1,
            frictionCoefficient: float=0,
            offset: float=0,
            glued: bool=False,
            relax = 1,
        ):
        super().__init__(value=value)

        self.type="implicitGapContact"

        self.patchType = patchType
        self.penaltyFactor = penaltyFactor
        self.penaltyFactorFriction = penaltyFactorFriction
        self.frictionCoefficient = frictionCoefficient
        self.relax = relax
        self.offset = offset
        self.glued = glued

    @property
    def patchType(self):
        return self._patchType

    @patchType.setter
    def patchType(self, patchType) -> None:
        check_type("patchType", patchType, str)
        check_value("patchType", patchType, _REGION_COUPLED_TYPES)
        self._patchType = patchType
        self.__setitem__('patchType', self.patchType)


    @property
    def penaltyFactor(self):
        return self._penaltyFactor

    @penaltyFactor.setter
    def penaltyFactor(self, penaltyFactor) -> None:
        check_type("penaltyFactor", penaltyFactor, (float, int))
        self._penaltyFactor = penaltyFactor
        self.__setitem__('penaltyFactor', self.penaltyFactor)


    @property
    def penaltyFactorFriction(self):
        return self._penaltyFactorFriction

    @penaltyFactorFriction.setter
    def penaltyFactorFriction(self, penaltyFactorFriction) -> None:
        check_type("penaltyFactorFriction", penaltyFactorFriction, (float, int))
        self._penaltyFactorFriction = penaltyFactorFriction
        self.__setitem__('penaltyFactorFriction', self.penaltyFactorFriction)


    @property
    def frictionCoefficient(self):
        return self._frictionCoefficient

    @frictionCoefficient.setter
    def frictionCoefficient(self, frictionCoefficient) -> None:
        check_type("frictionCoefficient", frictionCoefficient, (float, int))
        self._frictionCoefficient = frictionCoefficient
        self.__setitem__('frictionCoefficient', self.frictionCoefficient)


    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (float, int))
        self._relax = relax
        self.__setitem__('relax', self.relax)


    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, offset) -> None:
        check_type("offset", offset, (float, int))
        self._offset = offset
        self.__setitem__('offset', self.offset)


    @property
    def glued(self):
        return self._glued

    @glued.setter
    def glued(self, glued) -> None:
        check_type("glued", glued, bool)
        self._glued = glued
        self.__setitem__('glued', self.glued)