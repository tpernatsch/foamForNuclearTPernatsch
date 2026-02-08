from foamForNuclear.boundaryConditions.boundaryCondition import Patch
from foamForNuclear.boundaryConditions.tractionDisplacement import TractionDisplacement
from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import Table, Vector, List
from foamForNuclear.timeProfile import OffbeatTimeProfile


class PlenumSpringPressure(TractionDisplacement):
    """
    Coupled boundary condition that applies a **plenum spring pressure** to the
    displacement field.

    The `plenumSpringPressure` fvPatchField is selected in the patch subdictionary
    inside the `boundaryField` of the displacement field. It is intended to represent
    a mechanical spring in the plenum region acting between the top of the fuel column
    and the top cap / cladding top surface.

    The spring force is defined by a spring modulus and a pre-compression. The spring
    acts between two groups of patches: `fuelTopPatches` (bottom end of the spring)
    and `topCapInnerPatches` (top end of the spring). Optional external traction and/or
    pressure can also be applied through the inherited `tractionDisplacement` interface.

    Options
    -------
    value : list | Vector
        Initial displacement value on the patch.
        (required: True)

    springModulus : scalar
        Spring modulus in N/m, defining the stiffness of the plenum spring.
        (default: 3500.0; required: False)

    springPreCompression : scalar
        Pre-loading (pre-compression) of the spring in meters.
        (default: 0.0; required: False)

    fuelTopPatches : list[word]
        Patches attached to the bottom end of the spring (typically the top surface
        of the fuel column).
        (default: [fuelTop]; required: False)

    topCapInnerPatches : list[word]
        Patches attached to the top end of the spring (inner top-cap surface if present,
        or the top annular cladding surface in an open-cladding model).
        (default: [cladTop]; required: False)

    tractionList : OffbeatTimeProfile
        Time-dependent traction table (Pa). Used by the underlying `tractionDisplacement`
        interface.
        (required: False)

    traction : Vector
        Fixed traction vector (Pa). Used only if `tractionList` is absent.
        (required: False)

    pressureList : OffbeatTimeProfile
        Time-dependent pressure table (Pa). Used by the underlying `tractionDisplacement`
        interface.
        (required: False)

    pressure : scalar
        Fixed pressure (Pa). Used only if `pressureList` is absent.
        (required: False)

    planeStrain : bool
        Activate the plane strain approximation for the normal stress at the boundary.
        Cannot be used together with `flatSurface`.
        (default: False; required: False)

    flatSurface : bool
        Activate the flat surface approximation for the normal stress at the boundary.
        Cannot be used together with `planeStrain`.
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
        Relaxation factor for gradient updates.
        (default: 1; required: False)
    """

    def __init__(
            self,
            value,
            springModulus: float | int=3.5e3,
            springPreCompression: float | int=0.0,
            fuelTopPatches: list[str]=["fuelTop"],
            topCapInnerPatches: list[str]=["cladTop"],
            tractionList : OffbeatTimeProfile=None,
            traction: Vector=None,
            pressureList : OffbeatTimeProfile=None,
            pressure: float=None,
            planeStrain: bool=False,
            flatSurface: bool=False,
            fixedSpring: bool=False,
            fixedSpringModulus: float=None,
            dashpotModulus: float=None,
            relax: float=1
        ):
        super().__init__(
            value=value,
            tractionList=tractionList,
            traction=traction,
            pressureList=pressureList,
            pressure=pressure,
            planeStrain=planeStrain,
            flatSurface=flatSurface,
            fixedSpring=fixedSpring,
            fixedSpringModulus=fixedSpringModulus,
            dashpotModulus=dashpotModulus,
            relax=relax,
            )

        self.type = "plenumSpringPressure"
        self.springModulus = springModulus
        self.springPreCompression = springPreCompression
        self.fuelTopPatches = fuelTopPatches
        self.topCapInnerPatches = topCapInnerPatches

    @property
    def springModulus(self):
        return self._springModulus

    @springModulus.setter
    def springModulus(self, springModulus) -> None:
        check_type("springModulus", springModulus, (int, float))
        self._springModulus = springModulus
        self.__setitem__('springModulus', springModulus)

    @property
    def springPreCompression(self):
        return self._springPreCompression

    @springPreCompression.setter
    def springPreCompression(self, springPreCompression) -> None:
        check_type("springPreCompression", springPreCompression, (int, float))
        self._springPreCompression = springPreCompression
        self.__setitem__('springPreCompression', springPreCompression)

    @property
    def fuelTopPatches(self):
        return self._fuelTopPatches

    @fuelTopPatches.setter
    def fuelTopPatches(self, fuelTopPatches) -> None:
        check_type("fuelTopPatches", fuelTopPatches, list)
        self._fuelTopPatches = List(fuelTopPatches)
        self.__setitem__('fuelTopPatches', self.fuelTopPatches)

    @property
    def topCapInnerPatches(self):
        return self._topCapInnerPatches

    @topCapInnerPatches.setter
    def topCapInnerPatches(self, topCapInnerPatches) -> None:
        check_type("topCapInnerPatches", topCapInnerPatches, list)
        self._topCapInnerPatches = List(topCapInnerPatches)
        self.__setitem__('topCapInnerPatches', self.topCapInnerPatches)