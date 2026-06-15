from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Table, Vector
from foamForNuclear.timeProfile import OffbeatTimeProfile


class UnilateralContact(Patch):
    """
    Unilateral contact boundary condition for displacement, enforcing **non-penetration**
    against a contact plane.

    The `unilateralContact` fvPatchField is selected in the patch subdictionary
    inside the `boundaryField` of the displacement field. The contact constraint is
    enforced using a penalty formulation acting only in compression (no tensile
    resistance).

    A linear ramp can be applied to the contact stiffness as a function of the
    penetration depth to improve robustness and convergence.

    Options
    -------
    value : list | Vector
        Initial displacement value on the patch.
        (required: True)

    traction : Vector
        Fixed traction vector applied on the patch (Pa).
        (required: False)

    pressure : scalar
        Fixed pressure applied on the patch (Pa).
        (required: False)

    contactRelRampWidth : scalar
        Relative (with respect to the local cell size) width of the linear ramp
        used to progressively activate the contact stiffness.
        If not specified, the code interprets this as a small value.
        (required: False)

    contactRelOffset : scalar
        Relative (with respect to the local cell size) penetration tolerance used
        to define the contact plane.
        If not specified, the code interprets this as a small value.
        (required: False)

    penaltyScaleFactor : scalar
        Penalty scale factor controlling the contact stiffness.
        (default: 1.0; required: False)

    planeNormal : list | Vector
        Normal vector to the contact plane pointing OUTWARD from the solid
        (the direction of potential opening / contact approach). If not
        specified, the outward patch normal is used.
        (required: False)

    planePoint : list | Vector
        A point lying on the contact plane. If not specified, the code determines
        it automatically from the patch geometry.
        (required: False)

    relax : scalar
        Relaxation factor for stiffness weight updates.
        (default: 1.0; required: False)
    """

    def __init__(
            self,
            value,
            traction: Vector=None,
            pressure: float=None,
            contactRelRampWidth: float=None,
            contactRelOffset: float=None,
            penaltyScaleFactor: float=1.0,
            planeNormal: list|Vector=None,
            planePoint: list|Vector=None,
            relax: float=1.0,
        ):
        super().__init__(type="unilateralContact", value=value)

        self.traction = traction
        self.pressure = pressure
        self.contactRelRampWidth = contactRelRampWidth
        self.contactRelOffset = contactRelOffset
        self.penaltyScaleFactor = penaltyScaleFactor
        self.planeNormal = planeNormal
        self.planePoint = planePoint
        self.relax = relax

    @property
    def traction(self):
        return self._traction

    @traction.setter
    def traction(self, traction) -> None:
        check_type("traction", traction, Vector, none_ok=True)
        self._traction = traction
        if (traction is not None):
            self.__setitem__('traction', f"uniform {traction}")

    @property
    def pressure(self):
        return self._pressure

    @pressure.setter
    def pressure(self, pressure) -> None:
        check_type("pressure", pressure, (int, float), none_ok=True)
        self._pressure = pressure
        if (pressure is not None):
            self.__setitem__('pressure', f"uniform {pressure}")


    @property
    def contactRelRampWidth(self):
        return self._contactRelRampWidth

    @contactRelRampWidth.setter
    def contactRelRampWidth(self, contactRelRampWidth) -> None:
        check_type("contactRelRampWidth", contactRelRampWidth, (int, float), none_ok=True)
        self._contactRelRampWidth = contactRelRampWidth
        if contactRelRampWidth is not None:
            self.__setitem__('contactRelRampWidth', contactRelRampWidth)


    @property
    def contactRelOffset(self):
        return self._contactRelOffset

    @contactRelOffset.setter
    def contactRelOffset(self, contactRelOffset) -> None:
        check_type("contactRelOffset", contactRelOffset, (int, float), none_ok=True)
        self._contactRelOffset = contactRelOffset
        if contactRelOffset is not None:
            self.__setitem__('contactRelOffset', contactRelOffset)


    @property
    def penaltyScaleFactor(self):
        return self._penaltyScaleFactor

    @penaltyScaleFactor.setter
    def penaltyScaleFactor(self, penaltyScaleFactor) -> None:
        check_type("penaltyScaleFactor", penaltyScaleFactor, (int, float), none_ok=True)
        self._penaltyScaleFactor = penaltyScaleFactor
        if penaltyScaleFactor is not None:
            self.__setitem__('penaltyScaleFactor', penaltyScaleFactor)


    @property
    def planePoint(self):
        return self._planePoint

    @planePoint.setter
    def planePoint(self, planePoint) -> None:
        check_type("planePoint", planePoint, (list, Vector), none_ok=True)
        if planePoint is not None:
            self._planePoint = Vector(planePoint[0], planePoint[1], planePoint[2])
            self.__setitem__('planePoint', Vector(planePoint[0], planePoint[1], planePoint[2]))


    @property
    def planeNormal(self):
        return self._planeNormal

    @planeNormal.setter
    def planeNormal(self, planeNormal) -> None:
        check_type("planeNormal", planeNormal, (list, Vector), none_ok=True)
        if planeNormal is not None:
            self._planeNormal = Vector(planeNormal[0], planeNormal[1], planeNormal[2])
            self.__setitem__('planeNormal', Vector(planeNormal[0], planeNormal[1], planeNormal[2]))


    @property
    def relax(self):
        return self._relax

    @relax.setter
    def relax(self, relax) -> None:
        check_type("relax", relax, (int, float))
        self._relax = relax
        self.__setitem__('relax', relax)
