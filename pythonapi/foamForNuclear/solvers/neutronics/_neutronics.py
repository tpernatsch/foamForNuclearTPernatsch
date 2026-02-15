from __future__ import annotations

from copy import copy
import os

import attrs as attr
from attrs import define, field, validators as v
from foamForNuclear._attrs_tools import auto_type_validator

from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.externalSource import ExternalSource
from foamForNuclear.numerics import fvSchemes
from foamForNuclear.numerics import fvSolution, fvSolutionSolver
from foamForNuclear.mesh import Mesh
from foamForNuclear.nuclearData import ControlRodMove, NuclearData, NuclearDataState
from foamForNuclear.quadratureSet import QuadratureSet
from .._solvers import Solver
from foamForNuclear.timeFolder import TimeFolder


_NEUTRONICS_SOLVER_TYPES = {
    "diffusionNeutronics", "SP3Neutronics", "SNNeutronics", "adjointDiffusion",
    "pointKinetics"
}

class NeutronTransportOptions(OpenFOAMDict):
    """
    Neutron transport options. Handle acceleration flags and internal neutronics
    convergence criteria.

    Parameters
    ----------
    integralPredictor : bool
        Integral neutron balance made at each time step to predict fluxes at
        next step (can be unstable), (default `False`)
    implicitPredictor : bool
        Implicit soluion of all energy groups with explicit diffusion term
        (default `False`)
    ROMAcceleration : bool
        Reduced Order Method flag (default `False`)
    aitkenAcceleration : bool
        Aitken acceleration flag (default `False`)
    neutronIterationResidual : float
        Required accuracy for the coupling of different energy groups
        (default `1e-6`)
    maxNeutronIterations : int
        Up to 3-400 if no acceleration techniques (default `300`)
    """

    def __init__(
            self,
            integralPredictor: bool=False,
            implicitPredictor: bool=False,
            ROMAcceleration: bool=False,
            aitkenAcceleration: bool=False,
            neutronIterationResidual: float=1e-6,
            maxNeutronIterations: int=300,
        ):
        super().__init__()
        self.integralPredictor = integralPredictor
        self.implicitPredictor = implicitPredictor
        self.ROMAcceleration = ROMAcceleration
        self.aitkenAcceleration = aitkenAcceleration
        self.neutronIterationResidual = neutronIterationResidual
        self.maxNeutronIterations = maxNeutronIterations

    @property
    def integralPredictor(self):
        return self._integralPredictor

    @integralPredictor.setter
    def integralPredictor(self, integralPredictor) -> None:
        check_type("integralPredictor", integralPredictor, bool)
        self._integralPredictor = integralPredictor
        self.__setitem__('integralPredictor', integralPredictor)

    @property
    def implicitPredictor(self):
        return self._implicitPredictor

    @implicitPredictor.setter
    def implicitPredictor(self, implicitPredictor) -> None:
        check_type("implicitPredictor", implicitPredictor, bool)
        self._implicitPredictor = implicitPredictor
        self.__setitem__('implicitPredictor', implicitPredictor)

    @property
    def ROMAcceleration(self):
        return self._ROMAcceleration

    @ROMAcceleration.setter
    def ROMAcceleration(self, ROMAcceleration) -> None:
        check_type("ROMAcceleration", ROMAcceleration, bool)
        self._ROMAcceleration = ROMAcceleration
        self.__setitem__('ROMAcceleration', ROMAcceleration)

    @property
    def aitkenAcceleration(self):
        return self._aitkenAcceleration

    @aitkenAcceleration.setter
    def aitkenAcceleration(self, aitkenAcceleration) -> None:
        check_type("aitkenAcceleration", aitkenAcceleration, bool)
        self._aitkenAcceleration = aitkenAcceleration
        self.__setitem__('aitkenAcceleration', aitkenAcceleration)

    @property
    def neutronIterationResidual(self):
        return self._neutronIterationResidual

    @neutronIterationResidual.setter
    def neutronIterationResidual(self, neutronIterationResidual) -> None:
        check_type("neutronIterationResidual", neutronIterationResidual, float)
        check_positive("neutronIterationResidual", neutronIterationResidual)
        self._neutronIterationResidual = neutronIterationResidual
        self.__setitem__('neutronIterationResidual', neutronIterationResidual)

    @property
    def maxNeutronIterations(self):
        return self._maxNeutronIterations

    @maxNeutronIterations.setter
    def maxNeutronIterations(self, maxNeutronIterations) -> None:
        check_type("maxNeutronIterations", maxNeutronIterations, int)
        check_positive("maxNeutronIterations", maxNeutronIterations)
        self._maxNeutronIterations = maxNeutronIterations
        self.__setitem__('maxNeutronIterations', maxNeutronIterations)


@define(
    slots=True,
    on_setattr=[attr.setters.convert, attr.setters.validate],
    field_transformer=auto_type_validator,
    repr=False
)
class NeutronicsSolver(Solver):
    """
    Neutronics solver parameter object.

    Parameters
    ----------
    region : str
        Name of the neutronics region.
    solver : str {"diffusionNeutronics", "SP3Neutronics", "SNNeutronics", "adjointDiffusion", "pointKinetics"}
        Neutronics solver name.
    eigenvalueNeutronics : bool
        Eigenvalue calculation flag (default `True`).
    externalSourceNeutronics : bool
        External source calculation flag (default `False`).
    keff : float
        Initial eigenvalue (default `1`).
    power : float
        Initial power or normalization power (default `1` W).
    subcriticalIndex : float
        Optional: Subcritical index to be used only in subcritical
        point-kinetics (default `None`). If defined, don't define `ksrc`
    ksrc : float
        Optional: ksrc to be used only in subcritical point-kinetics
        (default `None`). If defined, don't define `subcriticalIndex`.
    timeFolder : TimeFolder
        Object representing the time folder (default `None`).
    mesh : Mesh
        Mesh object (default `None`).
    isMeshDeformation : bool
        Flag to allow mesh deformation based on a displacement field
        (default `False`).
    displacementFieldName : str
        Name of the displacement field (default `disp`).
    isSetFvSolutionToDefault : bool
        Flag to set fvSolution with default parameters (default `True`).
    isSetFvSchemesToDefault : bool
        Flag to set fvSchemes with default parameters (default `True`).

    Attributes
    ----------
    eigenvalueNeutronics : bool
        Eigenvalue calculation flag (default `True`).
    externalSourceNeutronics : bool
        External source calculation flag (default `False`).
    polyharmonicSplineMode : int
        Polyharmonic spline function mode (1. `r`; 2. `r^2 ln(r)`; 3. `r^3`;
        4. `r^4 ln(r)`) (default `1`).
    ScNo : float
        Schmidt number for diffusion of precursors (default `1`).
    axialOrientation : Vector
        Axial orientation of the axial expansion (default `None`).
    adjustDiscFactors : bool
        Flag to apply the discontinuity factor adjustement (default `False`).
    useGivenDiscFactors : bool
        Flag to use homogeneous discontinuity factors provided in state/zone
        (default `False`).
    legendreMoments : int
        Number of Legendre moments, used in SN solver (default `None`).
    isLowMemory : bool
        Flag for low memory foot print, useful for SN solver (defaulf `None`).
    fastNeutrons : bool
        Flag to change Doppler feedback correlation in the point-kinetics
        (default `True`).
    liquidFuel : bool
        Flag to specify that the fuel is liauid (default `False`).
    keff : float
        Initial eigenvalue (default `1`).
    power : float
        Initial power or normalization power (default `1` W).
    neutronTransportOptions : NeutronTransportOptions
    nuclearData : NuclearData
    quadratureSet : QuadratureSet
    externalSource : ExternalSource
    controlRodMove : ControlRodMove
    """
    region: str = ""
    solver: str = field(
        default="none",
        validator=v.and_(v.instance_of(str), v.in_(_NEUTRONICS_SOLVER_TYPES))
    )
    eigenvalueNeutronics: bool = True
    externalSourceNeutronics: bool = False
    polyharmonicSplineMode: int = 1
    ScNo: float | int = 1
    axialOrientation: Vector | None = None
    adjustDiscFactors: bool = False
    useGivenDiscFactors: bool = False
    legendreMoments: int | None = None
    isLowMemory: bool | None = None
    fastNeutrons: bool | None = None
    liquidFuel: bool = False
    keff: float | int = field(
        default=1,
        validator=v.and_(v.instance_of((float, int)), v.ge(0.0))
    )
    power: float | int = field(
        default=1,
        validator=v.and_(v.instance_of((float, int)), v.ge(0.0))
    )
    subcriticalIndex: float | int | None = None
    ksrc: float | int | None = None
    removeBaffles: bool = False
    timeFolder: TimeFolder | None = None
    mesh: Mesh | None = None
    isMeshDeformation: bool = False
    displacementFieldName: str = "disp"
    isSetFvSolutionToDefault: bool = True
    isSetFvSchemesToDefault: bool = True


    neutronTransportOptions: NeutronTransportOptions = field(factory=NeutronTransportOptions)
    nuclearData: NuclearData = field(init=False)
    quadratureSet: QuadratureSet = field(init=False)
    externalSource: ExternalSource = field(init=False)
    controlRodMove: ControlRodMove = field(init=False)


    def __attrs_post_init__(self):
        super().__attrs_post_init__()

        self.nuclearData = NuclearData(region=self.region)
        self.quadratureSet = QuadratureSet(region=self.region)
        self.externalSource = ExternalSource(region=self.region)
        self.controlRodMove = ControlRodMove(region=self.region)

        if self.isSetFvSolutionToDefault:
            self.set_fvSolution_default()
        if self.isSetFvSchemesToDefault:
            self.set_fvSchemes_default()


    def add_state(self, state: NuclearDataState):
        """
        Shortcut to add a nuclear data state in nuclearData.
        """
        self.nuclearData.add_state(state)


    def set_fvSolution_default(self):
        self.fvSolution = fvSolution()

        self.fvSolution.solvers["\"prec.*\""] = fvSolutionSolver(
            solver="PBiCG",
            preconditioner="DILU",
            tolerance=1e-6,
            relTol=1e-3,
        )
        self.fvSolution.solvers["\"precStar.*\""] = copy(self.fvSolution.solvers["\"prec.*\""])
        self.fvSolution.solvers["\"adjoint_prec.*\""] = copy(self.fvSolution.solvers["\"prec.*\""])

        self.fvSolution.solvers["\"flux.*\""] = fvSolutionSolver(
            solver="PCG",
            preconditioner="DIC",
            tolerance=1e-6,
            relTol=1e-3,
        )
        self.fvSolution.solvers["\"adjoint_flux.*\""] = copy(self.fvSolution.solvers["\"flux.*\""])

        self.fvSolution.solvers["\"angularFlux.*\""] = fvSolutionSolver(
            solver="PBiCGStab",
            preconditioner="DILU",
            tolerance=1e-7,
            relTol=1e-1,
        )


    def set_fvSchemes_default(self):
        self.fvSchemes = fvSchemes()

        self.fvSchemes.ddtSchemes['default'] = 'Euler'

        self.fvSchemes.gradSchemes['default'] = 'Gauss linear'

        self.fvSchemes.divSchemes['default'] = 'Gauss linear'
        self.fvSchemes.divSchemes['div(facePhi_,angularFlux_)'] = 'Gauss upwind'

        self.fvSchemes.laplacianSchemes['default'] = 'Gauss linear corrected'

        self.fvSchemes.interpolationSchemes['default'] = 'linear'

        self.fvSchemes.snGradSchemes['default'] = 'corrected'

        self.fvSchemes.fluxRequired['default'] = 'false'


    def get_required_fields(self) -> list[str]:
        """
        Return a list of required fields depending on the solver selected.
        They are labeled "MUST_READ"
        """
        if (self.solver == "diffusionNeutronics"):
            return(["defaultFlux"])
        if (self.solver == "SP3Neutronics"):
            return(["defaultFlux", "defaultFlux2"])
        if (self.solver == "SNNeutronics"):
            return(["defaultFlux"])
        if (self.solver == "adjointDiffusion"):
            return(["adjointDefaultPrec"])
        if (self.solver == "pointKinetics"):
            return([])


    def get_default_fields(self) -> list[str]:
        """
        Return a list of default fields that can be used as input depending on
        the solver selected. They are labeled "READ_IF_PRESENT"
        """
        common = [
            "powerDensity", "secondaryPowerDensity",
            "TFuel", "TClad", "TCool", "TStruct", "TStructMech"
        ]
        solverSpecific = []

        if (self.solver == "diffusionNeutronics"):
            solverSpecific = ["defaultExternalSourceFlux", "defaultPrec"]
        elif (self.solver == "SP3Neutronics"):
            solverSpecific = ["defaultPrec"]
        elif (self.solver == "SNNeutronics"):
            solverSpecific = ["defaultPrec"]
        elif (self.solver == "adjointDiffusion"):
            solverSpecific = []
        elif (self.solver == "pointKinetics"):
            solverSpecific = ["defaultPrec"]

        return(self.get_required_fields() + solverSpecific + common)


    def export_properties_to_openfoam(self):
        if (self.eigenvalueNeutronics and self.externalSourceNeutronics):
            msg = "eigenvalueNeutronics and externalSourceNeutronics are both true. Please select one of the two options"
            raise ValueError(msg)

        with open(f"constant/{self.region}/neutronicsProperties", 'w') as f:
            f.write(openfoamHeader)
            f.write(openfoamFileHeader("neutronicsProperties"))

            f.write(addParameter('eigenvalueNeutronics', self.eigenvalueNeutronics, isAddExtraLine=True))
            f.write(addParameter('externalSourceNeutronics', self.externalSourceNeutronics, isAddExtraLine=True))
            f.write(addParameter('polyharmonicSplineMode', self.polyharmonicSplineMode, isAddExtraLine=True))
            f.write(addParameter('ScNo', self.ScNo, isAddExtraLine=True))
            if (self.axialOrientation is not None):
                f.write(addParameter('axialOrientation', self.axialOrientation, isAddExtraLine=True))

            f.write(addParameter('adjustDiscFactors', self.adjustDiscFactors, isAddExtraLine=True, none_ok=False))
            f.write(addParameter('useGivenDiscFactors', self.useGivenDiscFactors, isAddExtraLine=True, none_ok=False))

            if (self.legendreMoments is not None):
                f.write(addParameter('legendreMoments', self.legendreMoments, isAddExtraLine=True, none_ok=False))
            if (self.isLowMemory is not None):
                f.write(addParameter('isLowMemory', self.isLowMemory, isAddExtraLine=True, none_ok=False))

            if (self.fastNeutrons is not None):
                f.write(addParameter("fastNeutrons", self.fastNeutrons, isAddExtraLine=True))
            f.write(addParameter('liquidFuel', self.liquidFuel, isAddExtraLine=True, none_ok=False))

            f.write(openfoamFooterLine)


    def export_uniform_to_openfoam(self):
        if (self.timeFolder is None):
            return

        self.timeFolder.export_to_openfoam()

        if (not os.path.exists(f"{self.timeFolder.time}/uniform")):
            os.mkdir(f"{self.timeFolder.time}/uniform")

        with open(f"{self.timeFolder.time}/uniform/reactorState", 'w') as f:
            f.write(openfoamHeader)
            f.write(openfoamFileHeader("reactorState"))

            f.write(addParameter('keff', self.keff, isAddExtraLine=True))
            f.write(addParameter('power', self.power, isAddExtraLine=True))

            if (self.ksrc is not None):
                f.write(addParameter('ksrc', self.ksrc, isAddExtraLine=True))
            if (self.subcriticalIndex is not None):
                f.write(addParameter('subcriticalIndex', self.subcriticalIndex, isAddExtraLine=True))

            f.write(openfoamFooterLine)


    def export_to_openfoam(self):
        self.create_folders()

        for obj in [
            self.nuclearData,
            self.quadratureSet,
            self.externalSource,
            self.controlRodMove,
            self.fvSolution,
            self.fvSchemes
        ]:
            obj.region = self.region

        self.export_properties_to_openfoam()
        self.export_uniform_to_openfoam()
        self.nuclearData.export_to_openfoam()

        if (self.quadratureSet is not None and self.solver == 'SNNeutronics'):
            self.quadratureSet.export_to_openfoam()

        if (self.externalSourceNeutronics and not self.externalSource.is_empty):
            self.externalSource.export_to_openfoam()

        if (not self.controlRodMove.is_empty):
            self.controlRodMove.export_to_openfoam()

        self.fvSolution.extraDict['neutronTransport'] = self.neutronTransportOptions

        return super().export_to_openfoam()
