from copy import copy

from foamForNuclear.common import *
from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import OpenFOAMDict
from foamForNuclear.dragModels import DragModel
from foamForNuclear.fvSchemes import fvSchemes
from foamForNuclear.heatTransferModels import HeatTransferModel
from foamForNuclear.phaseProperties import PhaseProperties, StructureProperty
from foamForNuclear.regimeMapModels import RegimeMapModel
from foamForNuclear.thermophysicalProperty import BaseThermophysicalProperty
from foamForNuclear.fvSolution import fvSolution, fvSolutionSolver
from foamForNuclear.solver import Solver
from foamForNuclear.transportProperties import TransportProperties
from foamForNuclear.turbulenceProperties import TurbulenceProperties


_THERMAL_HYDRAULICS_SOLVER_TYPES = {
    "onePhase", "twoPhase", "compressibleInterFoam", "rhoPimpleFoam"
}
_MOMENTUM_MODE_TYPES = {
    "cellCentered", "cellCenteredFaceReconstruction", "faceCentered"
}
_CONTINUITY_ERROR_COMPENSATION_MODE_TYPES = {"none", "Su", "Sp", "SuSp"}



class ThermalHydraulicsSolver(Solver):
    """
    Thermal-hydraulics solver parameter object.

    Parameters
    ----------
    region : str
        Name of the thermal-hydraulics region
    solver : str
        Thermal-hydraulics solver name (`onePhase`, `twoPhase`, or `compressibleInterFoam`)
    timeFolder : TimeFolder
        Object representing the time folder
    mesh : Mesh
        Mesh object
    isMeshDeformation : bool
        Flag to allow mesh deformation based on a displacement field (default `False`)
    displacementFieldName : str
        Name of the displacement field (default `disp`)
    isSetFvSolutionToDefault : bool
        Flag to set fvSolution with default parameters (default `True`)
    isSetFvSchemesToDefault : bool
        Flag to set fvSchemes with default parameters (default `True`)
    """
    def __init__(
            self,
            region = "",
            solver = "none",
            removeBaffles = False,
            timeFolder = None,
            mesh = None,
            isMeshDeformation = False,
            displacementFieldName = "disp",
            isSetFvSolutionToDefault: bool=True,
            isSetFvSchemesToDefault: bool=True
        ):
        super().__init__(region, solver, removeBaffles, timeFolder, mesh, isMeshDeformation, displacementFieldName)

        self.pimpleOptions: PimpleOptions = PimpleOptions()
        self.thermophysicalProperties: BaseThermophysicalProperty = BaseThermophysicalProperty(region=self.region)
        self.turbulenceProperties: TurbulenceProperties = TurbulenceProperties(region=self.region)
        self.phaseProperties: PhaseProperties = PhaseProperties(region=self.region)

        self.thermophysicalPropertiesSecondFluid: BaseThermophysicalProperty = BaseThermophysicalProperty(region=self.region)
        self.turbulencePropertiesSecondFluid: TurbulenceProperties = TurbulenceProperties(region=self.region)

        self.g = Vector(0, 0, -9.81)

        if (isSetFvSolutionToDefault):
            self.set_fvSolution_default()
        if (isSetFvSchemesToDefault):
            self.set_fvSchemes_default()


    @property
    def solver(self):
        return self._solver

    @solver.setter
    def solver(self, solver) -> None:
        if solver is not None:
            check_type("solver", solver, str)
            check_value("solver", solver, _THERMAL_HYDRAULICS_SOLVER_TYPES)
            self._solver = solver
        else:
            self._solver = ''

    @property
    def pimpleOptions(self):
        return self._pimpleOptions

    @pimpleOptions.setter
    def pimpleOptions(self, pimpleOptions) -> None:
        check_type("pimpleOptions", pimpleOptions, PimpleOptions)
        self._pimpleOptions = pimpleOptions

    @property
    def thermophysicalProperties(self):
        return self._thermophysicalProperties

    @thermophysicalProperties.setter
    def thermophysicalProperties(self, thermophysicalProperties) -> None:
        check_type("thermophysicalProperties", thermophysicalProperties, BaseThermophysicalProperty)
        self._thermophysicalProperties = thermophysicalProperties
        self._thermophysicalProperties.region = self._region

    @property
    def turbulenceProperties(self):
        return self._turbulenceProperties

    @turbulenceProperties.setter
    def turbulenceProperties(self, turbulenceProperties) -> None:
        check_type("turbulenceProperties", turbulenceProperties, TurbulenceProperties)
        self._turbulenceProperties = turbulenceProperties
        self._turbulenceProperties.region = self._region

    @property
    def phaseProperties(self):
        return self._phaseProperties

    @phaseProperties.setter
    def phaseProperties(self, phaseProperties) -> None:
        check_type("phaseProperties", phaseProperties, PhaseProperties)
        self._phaseProperties = phaseProperties
        self._phaseProperties.region = self._region

    @property
    def thermophysicalPropertiesSecondFluid(self):
        return self._thermophysicalPropertiesSecondFluid

    @thermophysicalPropertiesSecondFluid.setter
    def thermophysicalPropertiesSecondFluid(self, thermophysicalPropertiesSecondFluid) -> None:
        check_type("thermophysicalPropertiesSecondFluid", thermophysicalPropertiesSecondFluid, BaseThermophysicalProperty)
        self._thermophysicalPropertiesSecondFluid = thermophysicalPropertiesSecondFluid
        self._thermophysicalPropertiesSecondFluid.region = self._region

    @property
    def turbulencePropertiesSecondFluid(self):
        return self._turbulencePropertiesSecondFluid

    @turbulencePropertiesSecondFluid.setter
    def turbulencePropertiesSecondFluid(self, turbulencePropertiesSecondFluid) -> None:
        check_type("turbulencePropertiesSecondFluid", turbulencePropertiesSecondFluid, TurbulenceProperties)
        self._turbulencePropertiesSecondFluid = turbulencePropertiesSecondFluid
        self._turbulencePropertiesSecondFluid.region = self._region

    @property
    def g(self):
        return self._g

    @g.setter
    def g(self, g) -> None:
        check_type("g", g, Vector)
        self._g = g

    def add_drag_model(self, dragModel: DragModel):
        """
        Append a new drag model in the phaseProperties dictionnary
        """
        check_type("dragModel", dragModel, DragModel)
        self.phaseProperties.dragModels.append(dragModel)

    def add_heat_transfer_model(self, heatTransferModel: HeatTransferModel):
        """
        Append a new heat transfer model in the phaseProperties dictionnary
        """
        check_type("heatTransferModel", heatTransferModel, HeatTransferModel)
        self.phaseProperties.heatTransferModels.append(heatTransferModel)

    def add_regime_map_model(self, regimeMapModel: RegimeMapModel):
        """
        Append a new regime map model in the phaseProperties dictionnary
        """
        check_type("regimeMapModel", regimeMapModel, RegimeMapModel)
        self.phaseProperties.regimeMapModels.append(regimeMapModel)


    def add_structure_property(self, structureProperty: StructureProperty):
        """
        Append a new structure property in the phaseProperties dictionnary
        """
        check_type("structureProperty", structureProperty, StructureProperty)
        self.phaseProperties.structureProperties.append(structureProperty)

    def export_g_to_openfoam(self):
        with open(f"constant/{self.region}/g", 'w') as f:
            f.write(openfoamHeader)
            f.write(openfoamFileHeader("g", "uniformDimensionedVectorField"))

            f.write(addParameter('dimensions', "[0 1 -2 0 0 0 0]", isAddExtraLine=True))
            f.write(addParameter('value', self.g, isAddExtraLine=True))

            f.write(openfoamFooterLine)


    def set_fvSolution_default(self):
        self.fvSolution = fvSolution()

        self.fvSolution.solvers["p_rgh"] = fvSolutionSolver(
            solver="GAMG",
            smoother="DIC",
            tolerance=1e-7,
            relTol=1e-3,
        )
        self.fvSolution.solvers["p_rghFinal"] = fvSolutionSolver(
            solver="GAMG",
            smoother="DIC",
            tolerance=1e-7,
            relTol=0,
        )

        self.fvSolution.solvers["h"] = fvSolutionSolver(
            solver="PBiCG",
            preconditioner="DILU",
            tolerance=1e-7,
            relTol=1e-3,
        )
        self.fvSolution.solvers["k"] = copy(self.fvSolution.solvers["h"])
        self.fvSolution.solvers["epsilon"] = copy(self.fvSolution.solvers["h"])

        self.fvSolution.solvers['"e.*"'] = fvSolutionSolver(
            solver="smoothSolver",
            smoother="symGaussSeidel",
            tolerance=1e-08,
            relTol=0,
            minIter=0
        )

        self.fvSolution.solvers["hFinal"] = fvSolutionSolver(
            solver="PBiCG",
            preconditioner="DILU",
            tolerance=1e-7,
            relTol=0,
        )
        self.fvSolution.solvers["kFinal"] = copy(self.fvSolution.solvers["hFinal"])
        self.fvSolution.solvers["epsilonFinal"] = copy(self.fvSolution.solvers["hFinal"])


    def set_fvSchemes_default(self):
        self.fvSchemes = fvSchemes()

        self.fvSchemes.ddtSchemes['default'] = 'Euler'

        self.fvSchemes.gradSchemes['default'] = 'Gauss linear'

        self.fvSchemes.divSchemes['div(phi,alpha)'] = 'Gauss vanLeer'
        self.fvSchemes.divSchemes['div(phir,alpha)'] = 'Gauss vanLeer'
        self.fvSchemes.divSchemes['div(phirb,alpha)'] = 'Gauss linear'
        for scheme in [
            'div(rhoPhi,U)', 'div(rhoPhi,T)', 'div(rhoPhi,K)', 'div(phi,p)',
            'div(phi,k)'
        ]:
            self.fvSchemes.divSchemes[scheme] = 'Gauss upwind'
        self.fvSchemes.divSchemes['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'

        self.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

        self.fvSchemes.interpolationSchemes['default'] = 'linear'

        self.fvSchemes.snGradSchemes['default'] = 'uncorrected'


    def get_required_fields(self) -> list[str]:
        """
        Return a list of required fields depending on the solver selected.
        They are labeled "MUST_READ"
        """
        solverFields = []
        if (self.solver == "onePhase"):
            solverFields += ["p_rgh", "U", "T"]
        if (self.solver == "twoPhase"):
            solverFields += [
                "p_rgh", "U.<phase1>", "U.<phase2>", "T.<phase1>", "T.<phase2>"
            ]
        if (self.solver == "compressibleInterFoam"):
            solverFields += ["p_rgh", "U", "T"]

        if (self.turbulenceProperties.simulationType == "RAS"):
            if (self.turbulenceProperties.RASoptions.RASModel == "porousKEpsilon"):
                solverFields += ["k", "epsilon", "alphat", "nut"]

        return(solverFields)


    def get_default_fields(self) -> list[str]:
        """
        Return a list of default fields that can be used as input depending on
        the solver selected. They are labeled "READ_IF_PRESENT"
        """
        common = [
            "powerDensityStructure", "powerDensityLiquid",
            "magU", "alpha", "alphaPhi", "alphaRhoPhi", "alphaRhoMagU",
            "heatFlux.structure", "alpha.passiveStructure", "dgdt", "contErr",
            "T.passiveStructure", "T.fuelAvForNeutronics", "T.cladAvForNeutronics",
        ]
        solverSpecific = []

        if (self.solver == "onePhase"):
            solverSpecific = []
        elif (self.solver == "twoPhase"):
            solverSpecific = []
        elif (self.solver == "compressibleInterFoam"):
            solverSpecific = []

        return(self.get_required_fields() + solverSpecific + common)


    def export_to_openfoam(self):
        self.create_folders()

        # Only for two solver
        if (self.phaseProperties.isTwoPhase()):
            self.thermophysicalProperties.ext = self.phaseProperties.phaseNames[0]
            self.turbulenceProperties.ext = self.phaseProperties.phaseNames[0]
            self.thermophysicalPropertiesSecondFluid.ext += self.phaseProperties.phaseNames[1]
            self.turbulencePropertiesSecondFluid.ext += self.phaseProperties.phaseNames[1]

            self.thermophysicalPropertiesSecondFluid.region = self.region
            self.turbulencePropertiesSecondFluid.region = self.region

            self.thermophysicalPropertiesSecondFluid.export_to_openfoam()
            self.turbulencePropertiesSecondFluid.export_to_openfoam()

        self.thermophysicalProperties.region = self.region
        self.turbulenceProperties.region = self.region
        self.phaseProperties.region = self.region

        self.export_g_to_openfoam()
        self.thermophysicalProperties.export_to_openfoam()
        self.turbulenceProperties.export_to_openfoam()
        self.phaseProperties.export_to_openfoam()

        self.fvSolution.extraDict['PIMPLE'] = self.pimpleOptions

        return super().export_to_openfoam()


class CompressibleInterFoam(ThermalHydraulicsSolver):
    """
    OpenFOAM's solver for two compressible, non-isothermal immiscible fluids
    using a VOF (volume of fluid) phase-fraction based interface capturing
    approach.
    """
    def __init__(
            self,
            region="",
            removeBaffles=False,
            timeFolder=None,
            mesh=None,
            isMeshDeformation=False,
            displacementFieldName="disp",
            isSetFvSolutionToDefault: bool=True,
            isSetFvSchemesToDefault: bool=True
        ):
        super().__init__(
            region,
            "compressibleInterFoam",
            removeBaffles,
            timeFolder,
            mesh,
            isMeshDeformation,
            displacementFieldName,
            isSetFvSolutionToDefault
        )

        if (isSetFvSchemesToDefault):
            self.set_fvSchemes_default()

        if (isSetFvSolutionToDefault):
            self.set_fvSolution_default()


    def set_fvSchemes_default(self):
        self.fvSchemes = fvSchemes()

        self.fvSchemes.ddtSchemes['default'] = 'Euler'

        self.fvSchemes.gradSchemes['default'] = 'Gauss linear'

        self.fvSchemes.divSchemes['div(phi_,alpha)'] = 'Gauss vanLeer'
        self.fvSchemes.divSchemes['div(phirb,alpha)'] = 'Gauss linear'
        for scheme in [
            'div(rhoPhi,U)', 'div(rhoPhi,T)', 'div(rhoPhi,K)', 'div(phi,p)',
            'div(phi,k)'
        ]:
            self.fvSchemes.divSchemes[scheme] = 'Gauss upwind'
        self.fvSchemes.divSchemes['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'

        self.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

        self.fvSchemes.interpolationSchemes['default'] = 'linear'

        self.fvSchemes.snGradSchemes['default'] = 'uncorrected'

        del self.fvSchemes.fluxRequired['default']
        self.fvSchemes.fluxRequired['alpha.water'] = None
        self.fvSchemes.fluxRequired['p_rgh'] = None


    def set_fvSolution_default(self):
        self.fvSolution = fvSolution()

        self.fvSolution.solvers['"alpha.water.*"'] = fvSolutionSolver(
            solver="smoothSolver",
            smoother="symGaussSeidel",
            tolerance=1e-8,
            relTol=0,
            MULESCorr=True,
            nLimiterIter=5,
            nAlphaCorr=2,
            nAlphaSubCycles=1,
            cAlpha=1
        )

        self.fvSolution.solvers['".*(rho|rhoFinal)"'] = fvSolutionSolver(
            solver='diagonal'
        )

        self.fvSolution.solvers['"pcorr.*"'] = fvSolutionSolver(
            solver='PCG',
            preconditioner=OpenFOAMDict({
                'preconditioner': 'GAMG',
                'tolerance': 1e-05,
                'relTol': 0,
                'smoother': 'DICGaussSeidel',
            }),
            tolerance=1e-5,
            relTol=0,
            maxIter=100
        )

        self.fvSolution.solvers['p_rgh'] = fvSolutionSolver(
            solver='GAMG',
            tolerance=1e-7,
            relTol=0.01,
            smoother='DIC'
        )

        self.fvSolution.solvers['p_rghFinal'] = fvSolutionSolver(
            solver='PCG',
            preconditioner=OpenFOAMDict({
                'preconditioner': 'GAMG',
                'tolerance': 1e-05,
                'relTol': 0,
                'nVcycles': 2,
                'smoother': 'DICGaussSeidel',
                'nPreSweeps': 2
            }),
            tolerance=1e-7,
            relTol=0,
            maxIter=20
        )

        self.fvSolution.solvers['U'] = fvSolutionSolver(
            solver='smoothSolver',
            smoother='symGaussSeidel',
            tolerance=1e-6,
            relTol=0,
        )

        self.fvSolution.solvers['"(T|k|B|nuTilda).*"'] = fvSolutionSolver(
            solver='smoothSolver',
            smoother='symGaussSeidel',
            tolerance=1e-8,
            relTol=0,
        )

        self.pimpleOptions.nOuterCorrectors = 1
        self.pimpleOptions.nCorrectors = 3
        self.pimpleOptions.nNonOrthogonalCorrectors = 1
        self.pimpleOptions.momentumPredictor = False
        self.pimpleOptions.minNOuterCorrectors = None
        self.pimpleOptions.solveEnergy = None
        self.pimpleOptions.solveFluidMechanics = None
        self.pimpleOptions.momentumMode = None
        self.pimpleOptions.correctUntilConvergence = None
        self.pimpleOptions.porousInterfaceSharpness = None
        self.pimpleOptions.minMagU = None


    def export_to_openfoam(self):
        self.create_folders()

        self.thermophysicalProperties.region = self.region
        self.turbulenceProperties.region = self.region
        # self.phaseProperties.region = self.region

        self.export_g_to_openfoam()
        self.thermophysicalProperties.export_to_openfoam()
        self.turbulenceProperties.export_to_openfoam()
        # self.phaseProperties.export_to_openfoam()

        self.fvSolution.extraDict['PIMPLE'] = self.pimpleOptions

        # return super().export_to_openfoam()
        return(Solver.export_to_openfoam(self))



class PimpleOptions(OpenFOAMDict):
    """
    PIMPLE options

    Parameters
    ----------
    nCorrectors : int
    nOuterCorrectors : int
    solveEnergy : bool
        Solve the energy equation flag
    solveFluidMechanics : bool
        Solve fluid mechanics flag
    nNonOrthogonalCorrectors : int
    minNOuterCorrectors : int
    momentumMode : str
    correctUntilConvergence : bool
    porousInterfaceSharpness : float
    minMagU : float
        This limits the minimum magnitude of U that can be used in the
        calculation of the drag factor. Defaults to 0. Only (possibly) useful in
        certain two-phase scenarios.
    pMin : float
        Minimum pressure in Pa (default `10000`).
    pRefCell : int
        Cell index to apply the reference pressure (default `0`).
    pRefValue : float
        Reference pressure value (default `100000`).
    """

    def __init__(
            self,
            nCorrectors: int=2,
            nOuterCorrectors: int=2,
            solveEnergy: bool=True,
            solveFluidMechanics: bool=True,
            nNonOrthogonalCorrectors: int=0,
            minNOuterCorrectors: int=1,
            momentumPredictor: bool=None,
            momentumMode: str="cellCentered",
            correctUntilConvergence: bool=False,
            porousInterfaceSharpness: float=0.0,
            minMagU: float=0.0,
            pMin: float=10000,
            pRefCell: int=0,
            pRefValue: float=100000,
            continuityErrorCompensationMode: str=None,
            continuityErrorScaleFactor: float=None,
            partialEliminationMode: str=None,
            massTransferSafetyFactor: float=None,
            enthalpyStabilizationMode: str=None,
            oscillationLimiterFraction: float=None,
            maxTInterfaceDdt: float=None,
            correctPhi: bool=None
        ):
        super().__init__()
        self.nCorrectors = nCorrectors
        self.nNonOrthogonalCorrectors = nNonOrthogonalCorrectors
        self.nOuterCorrectors = nOuterCorrectors
        self.minNOuterCorrectors = minNOuterCorrectors
        self.solveEnergy = solveEnergy
        self.solveFluidMechanics = solveFluidMechanics
        self.momentumPredictor = momentumPredictor
        self.momentumMode = momentumMode
        self.correctUntilConvergence = correctUntilConvergence
        self.porousInterfaceSharpness = porousInterfaceSharpness
        self.minMagU = minMagU
        self.pMin = pMin
        self.pRefCell = pRefCell
        self.pRefValue = pRefValue
        self.continuityErrorCompensationMode = continuityErrorCompensationMode
        self.continuityErrorScaleFactor = continuityErrorScaleFactor
        self.partialEliminationMode = partialEliminationMode
        self.massTransferSafetyFactor = massTransferSafetyFactor
        self.enthalpyStabilizationMode = enthalpyStabilizationMode
        self.oscillationLimiterFraction = oscillationLimiterFraction
        self.maxTInterfaceDdt = maxTInterfaceDdt
        self.correctPhi = correctPhi

        self.residualControl = OpenFOAMListDict(OpenFOAMDict, name="residualControl")


    def __repr__(self, depth = 0):
        self.__setitem__('nCorrectors', self.nCorrectors)
        self.__setitem__('nNonOrthogonalCorrectors', self.nNonOrthogonalCorrectors)

        if (self.nOuterCorrectors is not None):
            self.__setitem__('nOuterCorrectors', self.nOuterCorrectors)
        if (self.solveEnergy is not None):
            self.__setitem__('solveEnergy', self.solveEnergy)
        if (self.solveFluidMechanics is not None):
            self.__setitem__('solveFluidMechanics', self.solveFluidMechanics)
        if (self.momentumMode is not None):
            self.__setitem__('momentumMode', self.momentumMode)
        if (self.correctUntilConvergence is not None):
            self.__setitem__('correctUntilConvergence', self.correctUntilConvergence)
        if (self.porousInterfaceSharpness is not None):
            self.__setitem__('porousInterfaceSharpness', self.porousInterfaceSharpness)
        if (self.minMagU is not None):
            self.__setitem__('minMagU', self.minMagU)
        if (self.pMin is not None):
            self.__setitem__('pMin', self.pMin)
        if (self.pRefCell is not None):
            self.__setitem__('pRefCell', self.pRefCell)
        if (self.pRefValue is not None):
            self.__setitem__('pRefValue', self.pRefValue)

        if (self.minNOuterCorrectors is not None):
            self.__setitem__('minNOuterCorrectors', self.minNOuterCorrectors)

        if (self.continuityErrorCompensationMode is not None):
            self.__setitem__('continuityErrorCompensationMode', self.continuityErrorCompensationMode)
        if (self.continuityErrorScaleFactor is not None):
            self.__setitem__('continuityErrorScaleFactor', self.continuityErrorScaleFactor)

        if (self.momentumPredictor is not None):
            self.__setitem__('momentumPredictor', self.momentumPredictor)
        if (self.partialEliminationMode is not None):
            self.__setitem__('partialEliminationMode', self.partialEliminationMode)
        if (self.massTransferSafetyFactor is not None):
            self.__setitem__('massTransferSafetyFactor', self.massTransferSafetyFactor)
        if (self.enthalpyStabilizationMode is not None):
            self.__setitem__('enthalpyStabilizationMode', self.enthalpyStabilizationMode)
        if (self.oscillationLimiterFraction is not None):
            self.__setitem__('oscillationLimiterFraction', self.oscillationLimiterFraction)
        if (self.maxTInterfaceDdt is not None):
            self.__setitem__('maxTInterfaceDdt', self.maxTInterfaceDdt)
        if (self.correctPhi is not None):
            self.__setitem__('correctPhi', self.correctPhi)

        if (not self.residualControl.is_empty):
            self.__setitem__('residualControl', self.residualControl)

        return super().__repr__(depth)



    @property
    def nCorrectors(self):
        return self._nCorrectors

    @nCorrectors.setter
    def nCorrectors(self, nCorrectors) -> None:
        check_type("nCorrectors", nCorrectors, int)
        check_positive("nCorrectors", nCorrectors)
        self._nCorrectors = nCorrectors

    @property
    def nNonOrthogonalCorrectors(self):
        return self._nNonOrthogonalCorrectors

    @nNonOrthogonalCorrectors.setter
    def nNonOrthogonalCorrectors(self, nNonOrthogonalCorrectors) -> None:
        check_type("nNonOrthogonalCorrectors", nNonOrthogonalCorrectors, int)
        check_positive("nNonOrthogonalCorrectors", nNonOrthogonalCorrectors)
        self._nNonOrthogonalCorrectors = nNonOrthogonalCorrectors

    @property
    def nOuterCorrectors(self):
        return self._nOuterCorrectors

    @nOuterCorrectors.setter
    def nOuterCorrectors(self, nOuterCorrectors) -> None:
        check_type("nOuterCorrectors", nOuterCorrectors, int, none_ok=True)
        if (nOuterCorrectors is not None):
            check_positive("nOuterCorrectors", nOuterCorrectors)
        self._nOuterCorrectors = nOuterCorrectors


    @property
    def minNOuterCorrectors(self):
        return self._minNOuterCorrectors

    @minNOuterCorrectors.setter
    def minNOuterCorrectors(self, minNOuterCorrectors) -> None:
        check_type("minNOuterCorrectors", minNOuterCorrectors, int, none_ok=True)
        if (minNOuterCorrectors is not None):
            check_positive("minNOuterCorrectors", minNOuterCorrectors)
        self._minNOuterCorrectors = minNOuterCorrectors

    @property
    def solveEnergy(self):
        return self._solveEnergy

    @solveEnergy.setter
    def solveEnergy(self, solveEnergy) -> None:
        check_type("solveEnergy", solveEnergy, bool, none_ok=True)
        self._solveEnergy = solveEnergy

    @property
    def solveFluidMechanics(self):
        return self._solveFluidMechanics

    @solveFluidMechanics.setter
    def solveFluidMechanics(self, solveFluidMechanics) -> None:
        check_type("solveFluidMechanics", solveFluidMechanics, bool, none_ok=True)
        self._solveFluidMechanics = solveFluidMechanics

    @property
    def momentumPredictor(self):
        return self._momentumPredictor

    @momentumPredictor.setter
    def momentumPredictor(self, momentumPredictor) -> None:
        check_type("momentumPredictor", momentumPredictor, bool, none_ok=True)
        self._momentumPredictor = momentumPredictor

    @property
    def momentumMode(self):
        return self._momentumMode

    @momentumMode.setter
    def momentumMode(self, momentumMode) -> None:
        check_type("momentumMode", momentumMode, str, none_ok=True)
        if momentumMode is not None:
            check_value("momentumMode", momentumMode, _MOMENTUM_MODE_TYPES)
        self._momentumMode = momentumMode

    @property
    def correctUntilConvergence(self):
        return self._correctUntilConvergence

    @correctUntilConvergence.setter
    def correctUntilConvergence(self, correctUntilConvergence) -> None:
        check_type("correctUntilConvergence", correctUntilConvergence, bool, none_ok=True)
        self._correctUntilConvergence = correctUntilConvergence

    @property
    def porousInterfaceSharpness(self):
        return self._porousInterfaceSharpness

    @porousInterfaceSharpness.setter
    def porousInterfaceSharpness(self, porousInterfaceSharpness) -> None:
        check_type("porousInterfaceSharpness", porousInterfaceSharpness, (int, float), none_ok=True)
        self._porousInterfaceSharpness = porousInterfaceSharpness

    @property
    def minMagU(self):
        return self._minMagU

    @minMagU.setter
    def minMagU(self, minMagU) -> None:
        check_type("minMagU", minMagU, (float, int), none_ok=True)
        self._minMagU = minMagU

    @property
    def pMin(self):
        return self._pMin

    @pMin.setter
    def pMin(self, pMin) -> None:
        check_type("pMin", pMin, (float, int))
        self._pMin = pMin

    @property
    def pRefCell(self):
        return self._pRefCell

    @pRefCell.setter
    def pRefCell(self, pRefCell) -> None:
        check_type("pRefCell", pRefCell, int)
        self._pRefCell = pRefCell

    @property
    def pRefValue(self):
        return self._pRefValue

    @pRefValue.setter
    def pRefValue(self, pRefValue) -> None:
        check_type("pRefValue", pRefValue, (float, int))
        self._pRefValue = pRefValue

    @property
    def continuityErrorCompensationMode(self):
        return self._continuityErrorCompensationMode

    @continuityErrorCompensationMode.setter
    def continuityErrorCompensationMode(self, continuityErrorCompensationMode) -> None:
        check_type("continuityErrorCompensationMode", continuityErrorCompensationMode, str, none_ok=True)

        if (continuityErrorCompensationMode is not None):
            check_value("continuityErrorCompensationMode", continuityErrorCompensationMode, _CONTINUITY_ERROR_COMPENSATION_MODE_TYPES)

        self._continuityErrorCompensationMode = continuityErrorCompensationMode

    @property
    def continuityErrorScaleFactor(self):
        return self._continuityErrorScaleFactor

    @continuityErrorScaleFactor.setter
    def continuityErrorScaleFactor(self, continuityErrorScaleFactor) -> None:
        check_type("continuityErrorScaleFactor", continuityErrorScaleFactor, (int, float), none_ok=True)
        self._continuityErrorScaleFactor = continuityErrorScaleFactor

    @property
    def partialEliminationMode(self):
        return self._partialEliminationMode

    @partialEliminationMode.setter
    def partialEliminationMode(self, partialEliminationMode) -> None:
        check_type("partialEliminationMode", partialEliminationMode, str, none_ok=True)
        self._partialEliminationMode = partialEliminationMode

    @property
    def massTransferSafetyFactor(self):
        return self._massTransferSafetyFactor

    @massTransferSafetyFactor.setter
    def massTransferSafetyFactor(self, massTransferSafetyFactor) -> None:
        check_type("massTransferSafetyFactor", massTransferSafetyFactor, (float, int), none_ok=True)
        self._massTransferSafetyFactor = massTransferSafetyFactor

    @property
    def enthalpyStabilizationMode(self):
        return self._enthalpyStabilizationMode

    @enthalpyStabilizationMode.setter
    def enthalpyStabilizationMode(self, enthalpyStabilizationMode) -> None:
        check_type("enthalpyStabilizationMode", enthalpyStabilizationMode, str, none_ok=True)
        self._enthalpyStabilizationMode = enthalpyStabilizationMode

    @property
    def oscillationLimiterFraction(self):
        return self._oscillationLimiterFraction

    @oscillationLimiterFraction.setter
    def oscillationLimiterFraction(self, oscillationLimiterFraction) -> None:
        check_type("oscillationLimiterFraction", oscillationLimiterFraction, (float, int), none_ok=True)
        self._oscillationLimiterFraction = oscillationLimiterFraction

    @property
    def maxTInterfaceDdt(self):
        return self._maxTInterfaceDdt

    @maxTInterfaceDdt.setter
    def maxTInterfaceDdt(self, maxTInterfaceDdt) -> None:
        check_type("maxTInterfaceDdt", maxTInterfaceDdt, (float, int), none_ok=True)
        self._maxTInterfaceDdt = maxTInterfaceDdt

    @property
    def correctPhi(self):
        return self._correctPhi

    @correctPhi.setter
    def correctPhi(self, correctPhi) -> None:
        check_type("correctPhi", correctPhi, bool, none_ok=True)
        self._correctPhi = correctPhi


    def add_residual_control_on_field(
            self,
            fieldName: str,
            tolerance: float,
            relTol: float,
            useFirstPISOInitialResidual: bool=None
        ):
        check_type("fieldName", fieldName, str)
        check_type("tolerance", tolerance, (float, int))
        check_type("relTol", relTol, (float, int))
        check_type("useFirstPISOInitialResidual", useFirstPISOInitialResidual, bool, none_ok=True)

        items = {
            "tolerance": tolerance,
            "relTol": relTol
        }
        if (useFirstPISOInitialResidual is not None):
            items["useFirstPISOInitialResidual"] = useFirstPISOInitialResidual

        self.residualControl.append(OpenFOAMDict(name=fieldName, items=items))
