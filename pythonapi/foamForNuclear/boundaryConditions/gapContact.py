from foamForNuclear.checkvalue import check_type, check_value
from .boundaryCondition import Patch
from .tractionDisplacement import TractionDisplacement


_REGION_COUPLED_TYPES = {"regionCoupledOFFBEAT"}


class GapContact(TractionDisplacement):
    """
    The `gapContact` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the displacement field.

    Parameters
    ----------
    patchType
        Specifies the patch field type. **Must be set to
        `regionCoupledOFFBEAT`.**
    penaltyFactor
        Penalty scale factor for the normal pressure. Controls the convergence
        and penetration balance. **Default: `0.1`.** Typical values: `0.01` to
        `1`.
    penaltyFactorFriction
        Penalty scale factor for the frictional forces. **Default: `0.1`.**
        Typical values: `0.01` to `1`.
    frictionCoefficient
        Coefficient of friction between surfaces. **Default: `0`.**
    relaxInterfacePressure
        Relaxation factor for the normal pressure update. **Default: `1`.**
        Typical values: `0.01` to `1`.
    relaxFriction
        Relaxation factor for the frictional force update. **Default: `1`.**
        Typical values: `0.01` to `1`.
    offset
        Offset added to the gap width calculation. A positive value increases
        the gap size. **Default: `0`.**
    rigidMasterNormal
        If `true`, the normal component is applied only on the slave side in
        rigid master mode. **Default: `false`.**
    rigidMasterFriction
        If `true`, the friction component is applied only on the slave side in
        rigid master mode. **Default: `false`.**

    Parameters in the patch subdictionary derived from `tractionDisplacement`:

    planeStrain
        Activates the plane strain approximation for the normal stress at the
        boundary. When enabled, the normal strain is assumed constant across the
        last layer of cells. **Default: `false`.**
    fixedSpring
        Activates a fixed spring-dashpot system for additional stability.
        **Default: `false`.**
    fixedSpringModulus
        Spring modulus in N/m. Required when `fixedSpring` is set to `true`.
    dashpotModulus
        Dashpot modulus in N/m. Required when `fixedSpring` is set to `true`.
    relax
        Relaxation factor for gradient updates. **Default: `1.0`.**
    value
        Initial displacement value (not stress).
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
