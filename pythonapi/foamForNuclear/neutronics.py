from copy import copy
import os

from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import *
from foamForNuclear.externalSource import ExternalSource
from foamForNuclear.fvSchemes import fvSchemes
from foamForNuclear.fvSolution import fvSolution, fvSolutionSolver
from foamForNuclear.mesh.mesh import Mesh
from foamForNuclear.nuclearData import ControlRodMove, NuclearData, NuclearDataState
from foamForNuclear.quadratureSet import QuadratureSet
from foamForNuclear.solver import Solver
from foamForNuclear.timeFolder import TimeFolder


_NEUTRONICS_SOLVER_TYPES = {
    "diffusionNeutronics", "SP3Neutronics", "SNNeutronics", "adjointDiffusion",
    "pointKinetics"
}


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
    def __init__(
            self,
            region: str="",
            solver: str="none",
            eigenvalueNeutronics: bool=True,
            externalSourceNeutronics: bool=False,
            keff: float=1,
            power: float=1,
            subcriticalIndex: float=None,
            ksrc: float=None,
            removeBaffles: bool=False,
            timeFolder: TimeFolder=None,
            mesh: Mesh=None,
            isMeshDeformation: bool=False,
            displacementFieldName: str="disp",
            isSetFvSolutionToDefault: bool=True,
            isSetFvSchemesToDefault: bool=True
        ):
        super().__init__(
            region, solver, removeBaffles, timeFolder, mesh, isMeshDeformation,
            displacementFieldName
        )

        self.eigenvalueNeutronics = eigenvalueNeutronics
        self.externalSourceNeutronics = externalSourceNeutronics
        self.keff = keff
        self.power = power

        self.subcriticalIndex = subcriticalIndex
        self.ksrc = ksrc

        self.neutronTransportOptions: NeutronTransportOptions = NeutronTransportOptions()
        self.nuclearData: NuclearData = NuclearData(region=region)
        self.quadratureSet: QuadratureSet = QuadratureSet(region=region)
        self.externalSource: ExternalSource = ExternalSource(region=region)
        self.controlRodMove: ControlRodMove = ControlRodMove(region=region)

        if (isSetFvSolutionToDefault):
            self.set_fvSolution_default()
        if (isSetFvSchemesToDefault):
            self.set_fvSchemes_default()


    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, solver) -> None:
        check_type("solver", solver, str)
        check_value("solver", solver, _NEUTRONICS_SOLVER_TYPES)
        self._solver = solver

    @property
    def eigenvalueNeutronics(self):
        return self._eigenvalueNeutronics

    @eigenvalueNeutronics.setter
    def eigenvalueNeutronics(self, eigenvalueNeutronics) -> None:
        check_type("eigenvalueNeutronics", eigenvalueNeutronics, bool)
        self._eigenvalueNeutronics = eigenvalueNeutronics

    @property
    def externalSourceNeutronics(self):
        return self._externalSourceNeutronics

    @externalSourceNeutronics.setter
    def externalSourceNeutronics(self, externalSourceNeutronics) -> None:
        check_type("externalSourceNeutronics", externalSourceNeutronics, bool)
        self._externalSourceNeutronics = externalSourceNeutronics

    @property
    def keff(self):
        return self._keff

    @keff.setter
    def keff(self, keff) -> None:
        check_type("keff", keff, (float, int))
        check_positive("keff", keff, is_strict=True)
        self._keff = keff

    @property
    def power(self):
        return self._power

    @power.setter
    def power(self, power) -> None:
        check_type("power", power, (float, int))
        check_positive("power", power)
        self._power = power

    @property
    def ksrc(self):
        return self._ksrc

    @ksrc.setter
    def ksrc(self, ksrc):
        check_type("ksrc", ksrc, (int, float), none_ok=True)
        self._ksrc = ksrc

    @property
    def subcriticalIndex(self):
        return self._subcriticalIndex

    @subcriticalIndex.setter
    def subcriticalIndex(self, subcriticalIndex):
        check_type("subcriticalIndex", subcriticalIndex, (int, float), none_ok=True)
        self._subcriticalIndex = subcriticalIndex

    @property
    def neutronTransportOptions(self):
        return self._neutronTransportOptions

    @neutronTransportOptions.setter
    def neutronTransportOptions(self, neutronTransportOptions) -> None:
        check_type("neutronTransportOptions", neutronTransportOptions, NeutronTransportOptions)
        self._neutronTransportOptions = neutronTransportOptions

    @property
    def nuclearData(self):
        return self._nuclearData

    @nuclearData.setter
    def nuclearData(self, nuclearData) -> None:
        check_type("nuclearData", nuclearData, NuclearData)
        self._nuclearData = nuclearData
        self._nuclearData.region = self.region

    @property
    def quadratureSet(self):
        return self._quadratureSet

    @quadratureSet.setter
    def quadratureSet(self, quadratureSet) -> None:
        check_type("quadratureSet", quadratureSet, QuadratureSet)
        self._quadratureSet = quadratureSet
        self._quadratureSet.region = self.region

    @property
    def externalSource(self):
        return self._externalSource

    @externalSource.setter
    def externalSource(self, externalSource) -> None:
        check_type("externalSource", externalSource, ExternalSource)
        self._externalSource = externalSource
        self._externalSource.region = self.region

    @property
    def controlRodMove(self):
        return self._controlRodMove

    @controlRodMove.setter
    def controlRodMove(self, controlRodMove) -> None:
        check_type("controlRodMove", controlRodMove, ControlRodMove)
        self._controlRodMove = controlRodMove
        self._controlRodMove.region = self.region


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

            f.write(openfoamFooterLine)


    def export_uniform_to_openfoam(self):
        if (self.timeFolder is not None):
            self.timeFolder.export_to_openfoam()

            if (not os.path.exists(f"{self.timeFolder.time}/uniform")):
                os.mkdir(f"{self.timeFolder.time}/uniform")

            with open(f"{self.timeFolder.time}/uniform/reactorState", 'w') as f:
                f.write(openfoamHeader)
                f.write(openfoamFileHeader("reactorState"))

                f.write(addParameter('keff', self.keff, isAddExtraLine=True))
                f.write(addParameter('pTarget', self.power, isAddExtraLine=True))

                if (self.ksrc is not None):
                    f.write(addParameter('ksrc', self.ksrc, isAddExtraLine=True))
                if (self.subcriticalIndex is not None):
                    f.write(addParameter('subcriticalIndex', self.subcriticalIndex, isAddExtraLine=True))

                f.write(openfoamFooterLine)


    def export_to_openfoam(self):
        self.create_folders()

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
