from foamForNuclear.checkvalue import check_type, check_value
from .boundaryCondition import Patch
from .tractionDisplacement import TractionDisplacement


_REGION_COUPLED_TYPES = {"regionCoupledOFFBEAT"}


class GapContact(TractionDisplacement):
    """
    Coupled contact boundary condition for displacement, based on a penalty method.

    The `gapContact` fvPatchField is selected in the patch subdictionary inside the
    `boundaryField` of the displacement field. It is intended to be used on a *pair*
    of `regionCoupledOFFBEAT` patches coupled via AMI.

    The normal contact pressure is enforced with a penalty approach. Optional friction
    can be activated via a Coulomb friction coefficient and a separate penalty scaling
    for tangential forces.

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
        between penetration and convergence.
        Typical values: 0.01 to 1.
        (default: 0.1; required: False)

    penaltyFactorFriction : scalar
        Penalty scale factor for the frictional (tangential) forces.
        Typical values: 0.01 to 1.
        (default: 0.1; required: False)

    frictionCoefficient : scalar
        Coulomb friction coefficient between surfaces.
        (default: 0; required: False)

    relaxInterfacePressure : scalar
        Relaxation factor for the normal pressure update.
        Typical values: 0.01 to 1.
        (default: 1; required: False)

    relaxFriction : scalar
        Relaxation factor for the frictional force update.
        Typical values: 0.01 to 1.
        (default: 1; required: False)

    offset : scalar
        Offset added to the computed gap width. A positive value increases the gap.
        (default: 0; required: False)

    rigidMasterNormal : bool
        If true, apply the normal component only on the slave side in rigid-master mode.
        (default: False; required: False)

    rigidMasterFriction : bool
        If true, apply the friction component only on the slave side in rigid-master mode.
        (default: False; required: False)

    planeStrain : bool
        Activate the plane strain approximation for the normal stress at the boundary.
        (default: False; required: False)

    fixedSpring : bool
        Enable a fixed spring-dashpot system for additional stability.
        (default: False; required: False)

    fixedSpringModulus : scalar
        Spring modulus in N/m (only if `fixedSpring` is true).
        (required: False)

    dashpotModulus : scalar
        Dashpot modulus in N/m (only if `fixedSpring` is true).
        (required: False)

    relax : scalar
        Relaxation factor for gradient updates. **It is not suggest to use this 
        relaxation factor. Prefer `relaxInterfacePressure` to improve convergence**.
        (default: 1; required: False)

    traction : Vector
        Initial traction value (from `tractionDisplacement` interface).
        (required: False)

    tractionList : list[Vector]
        Time-dependent traction specification (from `tractionDisplacement` interface).
        (required: False)

    pressure : scalar
        Initial pressure value (from `tractionDisplacement` interface).
        (required: False)

    pressureList : list[scalar]
        Time-dependent pressure specification (from `tractionDisplacement` interface).
        (required: False)
    """

    def __init__(
            self,
            value,
            patchType: str="regionCoupledOFFBEAT",
            penaltyFactor: float=0.1,
            penaltyFactorFriction: float=0.1,
            frictionCoefficient: float=0,
            relaxInterfacePressure: float=0,
            relaxFriction: float=1,
            offset: float=0,
            rigidMasterNormal: bool=False,
            rigidMasterFriction: bool=False,
            tractionList = None,
            traction = None,
            pressureList = None,
            pressure = None,
            fixedSpringModulus = None,
            dashpotModulus = None,
            planeStrain = False,
            flatSurface = False,
            fixedSpring = False,
            relax = 1
        ):
        super().__init__(
            value=value,
            tractionList=tractionList,
            traction=traction,
            pressureList=pressureList,
            pressure=pressure,
            fixedSpringModulus=fixedSpringModulus,
            dashpotModulus=dashpotModulus,
            planeStrain=planeStrain,
            flatSurface=flatSurface,
            fixedSpring=fixedSpring,
            relax=relax
        )

        self.type="gapContact"

        self.patchType = patchType
        self.penaltyFactor = penaltyFactor
        self.penaltyFactorFriction = penaltyFactorFriction
        self.frictionCoefficient = frictionCoefficient
        self.relaxInterfacePressure = relaxInterfacePressure
        self.relaxFriction = relaxFriction
        self.offset = offset
        self.rigidMasterNormal = rigidMasterNormal
        self.rigidMasterFriction = rigidMasterFriction


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
    def relaxInterfacePressure(self):
        return self._relaxInterfacePressure

    @relaxInterfacePressure.setter
    def relaxInterfacePressure(self, relaxInterfacePressure) -> None:
        check_type("relaxInterfacePressure", relaxInterfacePressure, (float, int))
        self._relaxInterfacePressure = relaxInterfacePressure
        self.__setitem__('relaxInterfacePressure', self.relaxInterfacePressure)


    @property
    def relaxFriction(self):
        return self._relaxFriction

    @relaxFriction.setter
    def relaxFriction(self, relaxFriction) -> None:
        check_type("relaxFriction", relaxFriction, (float, int))
        self._relaxFriction = relaxFriction
        self.__setitem__('relaxFriction', self.relaxFriction)


    @property
    def offset(self):
        return self._offset

    @offset.setter
    def offset(self, offset) -> None:
        check_type("offset", offset, (float, int))
        self._offset = offset
        self.__setitem__('offset', self.offset)


    @property
    def rigidMasterNormal(self):
        return self._rigidMasterNormal

    @rigidMasterNormal.setter
    def rigidMasterNormal(self, rigidMasterNormal) -> None:
        check_type("rigidMasterNormal", rigidMasterNormal, bool)
        self._rigidMasterNormal = rigidMasterNormal
        self.__setitem__('rigidMasterNormal', self.rigidMasterNormal)


    @property
    def rigidMasterFriction(self):
        return self._rigidMasterFriction

    @rigidMasterFriction.setter
    def rigidMasterFriction(self, rigidMasterFriction) -> None:
        check_type("rigidMasterFriction", rigidMasterFriction, bool)
        self._rigidMasterFriction = rigidMasterFriction
        self.__setitem__('rigidMasterFriction', self.rigidMasterFriction)
