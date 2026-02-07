from __future__ import annotations
from typing import Optional

from io import StringIO
import os

import attrs as attr
from attrs import define, field

import foamForNuclear
from foamForNuclear._attrs_tools import auto_type_validator, call_method_on_change, _propagate_region_to, ffn_define

from copy import copy

from foamForNuclear.common import *
from foamForNuclear.checkvalue import check_positive, check_type, check_value
from foamForNuclear.common import OpenFOAMDict
from foamForNuclear.porous_medium import pair_geometry
from foamForNuclear.porous_medium.drag import DragModel
from foamForNuclear.numerics import fvSchemes
from foamForNuclear.porous_medium.dispersed_diameter import DispersedDiameterModel
from foamForNuclear.porous_medium.phase_change import PhaseChangeModel
from foamForNuclear.porous_medium.two_phase_drag_multiplier import TwoPhaseDragMultiplierModel
from foamForNuclear.porous_medium.pair_geometry import PairGeometryModel, ContactPartitionModel, DispersionModel, InterfacialAreaDensityModel
from foamForNuclear.porous_medium.heat_transfer import HeatTransferModel
from foamForNuclear.porous_medium import Structure, Fluid
from foamForNuclear.porous_medium.regime_map import RegimeMapModel
from foamForNuclear.thermo import BaseThermophysicalProperty
from foamForNuclear.numerics import fvSolution, fvSolutionSolver
from .._solvers import Solver
from foamForNuclear.transport import TransportProperties
from foamForNuclear.turbulence import TurbulenceProperties



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
        check_type("momentumMode", momentumMode, str)
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


_THERMAL_HYDRAULICS_SOLVER_TYPES = {
    "onePhase", "twoPhase", "compressibleInterFoam", "rhoPimpleFoam"
}
_MOMENTUM_MODE_TYPES = {
    "cellCentered", "cellCenteredFaceReconstruction", "faceCentered"
}
_CONTINUITY_ERROR_COMPENSATION_MODE_TYPES = {"none", "Su", "Sp", "SuSp"}


@ffn_define
class FluidStructureModels(FoamForNuclearDict):
    """
    Models coupling the fluid with surrounding structures.
    """
    dragModels: list[DragModel] = field(factory=list)
    heatTransferModels: list[HeatTransferModel] = field(factory=list)
    pairGeometryModel: Optional[PairGeometryModel] = None


@ffn_define
class FluidFluidModels(FoamForNuclearDict):
    """
    Models for fluid–fluid interactions.
    """
    dragModels: list[DragModel] = field(factory=list)
    heatTransferModels: list[HeatTransferModel] = field(factory=list)
    twoPhaseDragMultiplierModel: Optional[TwoPhaseDragMultiplierModel] = None
    pairGeometryModel: Optional[PairGeometryModel] = None
    phaseChangeModel: Optional[PhaseChangeModel] = None


@define(
    slots=True,
    on_setattr=[attr.setters.validate, call_method_on_change("propagate_region", "region")],
    field_transformer=auto_type_validator,
    repr=False,
)
class ThermalHydraulicsSolver(Solver):
    """
    Base thermal-hydraulics solver parameters.

    Specialised solvers (one-phase, two-phase, compressibleInterFoam, ...)
    should inherit from this and only add their extra fields.
    """
    solver: str = field(default="none", init=False, on_setattr=attr.setters.frozen,)
    region: str = ""
    isSetFvSolutionToDefault: bool = True
    isSetFvSchemesToDefault: bool = True
    g: Vector = Vector(0.0, 0.0, -9.81)
    pimpleOptions: PimpleOptions = field(factory=PimpleOptions)
    pMin: int | float = 10000
    pRefCell: int  = 0
    pRefValue: int | float = 100000
    residualKd: int | float | None = None

    @solver.validator
    def _solver_allowed(self, _attr, value: str):
        if value not in _THERMAL_HYDRAULICS_SOLVER_TYPES:
            raise ValueError(f"Provided solver type \"{value}\" does not exists. Solver must be one of {sorted(_THERMAL_HYDRAULICS_SOLVER_TYPES)}")
        
    
    def propagate_region(self, new_value: str | None = None) -> None:
        super().propagate_region(new_value)
        if new_value is None:
            new_value = self.region
        _propagate_region_to(self, new_value, "thermophysicalProperties", "turbulenceProperties", 
                             "phaseProperties", "thermophysicalPropertiesSecondFluid", 
                             "turbulencePropertiesSecondFluid")

    def __attrs_post_init__(self):
        if self.isSetFvSolutionToDefault:
            self.set_fvSolution_default()
        if self.isSetFvSchemesToDefault:
            self.set_fvSchemes_default()

        super().__attrs_post_init__()
        _propagate_region_to(self, self.region, "thermophysicalProperties", "turbulenceProperties", 
                             "phaseProperties", "thermophysicalPropertiesSecondFluid", 
                             "turbulencePropertiesSecondFluid")


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

        # if (self.turbulenceProperties.simulationType == "RAS"):
        #     if (self.turbulenceProperties.RASoptions.RASModel == "porousKEpsilon"):
        #         solverFields += ["k", "epsilon", "alphat", "nut"]

        return(solverFields)


    def get_default_fields(self) -> list[str]:
        """
        Return a list of default fields that can be used as input depending on
        the solver selected. They are labeled "READ_IF_PRESENT"
        """
        common = [
            "powerDensityNeutronics", "powerDensityNeutronicsToLiquid",
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

    def export_g_to_openfoam(self):
        with open(f"constant/{self.region}/g", 'w') as f:
            f.write(openfoamHeader)
            f.write(openfoamFileHeader("g", "uniformDimensionedVectorField"))

            f.write(addParameter('dimensions', "[0 1 -2 0 0 0 0]", isAddExtraLine=True))
            f.write(addParameter('value', self.g, isAddExtraLine=True))

            f.write(openfoamFooterLine)

    def export_phase_properties_to_openfoam(self):
        filename = "phaseProperties"

        # --- pick folder / region path -----------------------------------------
        if self.region in ("", "region0"):
            out_dir = "constant"
        else:
            out_dir = f"constant/{self.region}"

        os.makedirs(out_dir, exist_ok=True)
        out_path = f"{out_dir}/{filename}"

        buf = StringIO()

        # header
        buf.write(openfoamHeader)
        buf.write(openfoamFileHeader(filename))

        # body: delegate to subclass
        self._write_phase_properties_body(buf)

        # footer
        buf.write(openfoamFooterLine)

        with open(out_path, "w", encoding="utf-8") as f:
            f.write(buf.getvalue())

    def _write_phase_properties_body(self, buf: StringIO):
        raise NotImplementedError

    def export_to_openfoam(self):
        return super().export_to_openfoam()


@define(
    slots=True,
    on_setattr=[attr.setters.validate, call_method_on_change("propagate_region", "region")],
    field_transformer=auto_type_validator,
    repr=False,
)
class OnePhaseThermalHydraulicsSolver(ThermalHydraulicsSolver):
    """
    One-phase thermal-hydraulics solver.

    High-level description:
      - one fluid
      - surrounding structures
      - fluid–structure interaction models
    """
    solver: str = field(default="onePhase", init=False, on_setattr=attr.setters.frozen,)
    fluid: Fluid = field(factory=Fluid)
    structures: list[Structure] = field(factory=list)
    fluid_structure: FluidStructureModels = field(factory=FluidStructureModels)
    regimeMapModels: list[RegimeMapModel] = field(factory=list)

    def _write_phase_properties_body(self, buf: StringIO):
        # Structures        
        structure_properties_dict = OpenFOAMListDict(
            name="structureProperties", expected_type=Structure, items=self.structures)
        buf.write(f"{structure_properties_dict!r}\n")

        # Regime maps
        regimes_dict = OpenFOAMListDict(
            name="regimeMapModels", expected_type=RegimeMapModel, items=self.regimeMapModels)
        buf.write(f"{regimes_dict!r}\n")

        # Physics models
        buf.write(f"physicsModels\n{{\n")

        # Drag
        fluid_structure_drag = OpenFOAMListDict(
            name="dragModels", expected_type=DragModel, items=self.fluid_structure.dragModels)
        buf.write(f"{fluid_structure_drag.__repr__(depth=1)}\n")

        # Heat transfer
        fluid_structure_ht = OpenFOAMListDict(
            name="heatTransferModels", expected_type=HeatTransferModel, items=self.fluid_structure.heatTransferModels)
        buf.write(f"{fluid_structure_ht.__repr__(depth=1)}\n")
        buf.write("}")

        # pMin, pRefCell, pRefValue, residualKd
        buf.write(f"pMin {self.pMin};\n")
        buf.write(f"pRefCell {self.pRefCell};\n")
        buf.write(f"pRefValue {self.pRefValue};\n")
        buf.write(f"residualKd {self.residualKd};\n")

    def export_to_openfoam(self):
        self.create_folders()

        self.fluid.turbulenceProperties.region = self.region
        self.fluid.turbulenceProperties.export_to_openfoam()

        self.fluid.thermophysicalProperties.region = self.region
        self.fluid.thermophysicalProperties.export_to_openfoam()

        self.export_g_to_openfoam()

        self.export_phase_properties_to_openfoam()

        self.fvSolution.extraDict['PIMPLE'] = self.pimpleOptions

        return super().export_to_openfoam()


@define(
    slots=True,
    on_setattr=[attr.setters.validate, attr.setters.convert, call_method_on_change("propagate_region", "region")],
    field_transformer=auto_type_validator,
    repr=False,
    kw_only=True,
)
class TwoPhaseThermalHydraulicsSolver(ThermalHydraulicsSolver):
    """
    Two-phase thermal-hydraulics solver.

    High-level description:
      - two fluids
      - surrounding structures
      - fluid(s)–structure interaction models
      - fluid–fluid interaction models
    """
    solver: str = field(default="twoPhase", init=False, on_setattr=attr.setters.frozen,)
    fluid1: Fluid = field(
        default=attr.Factory(lambda self: Fluid(name="fluid", stateOfMatter="liquid"), takes_self=True))
    fluid2: Fluid = field(
        default=attr.Factory(lambda self: Fluid(name="vapour", stateOfMatter="gas"), takes_self=True))
    structures: list[Structure] = field(factory=list)
    fluid1_structure: FluidStructureModels = field(factory=FluidStructureModels)
    fluid2_structure: FluidStructureModels = field(factory=FluidStructureModels)
    fluid_fluid: FluidFluidModels = field(factory=FluidFluidModels)
    regimeMapModels: list[RegimeMapModel] = field(factory=list)

    def _write_phase_properties_body(self, buf: StringIO):
        # Fluids & structures
        buf.write(f"fluid1 {self.fluid1.name};\n")
        buf.write(f"fluid2 {self.fluid2.name};\n\n")

        buf.write(f"{self.fluid1.name}Properties{self.fluid1._as_openfoam_dict().export_body_to_foam()}\n")
        buf.write(f"{self.fluid2.name}Properties{self.fluid2._as_openfoam_dict().export_body_to_foam()}\n")
        
        structure_properties_dict = OpenFOAMListDict(
            name="structureProperties", expected_type=Structure, items=self.structures)
        buf.write(f"{structure_properties_dict!r}\n")

        # Regime maps
        regimes_dict = OpenFOAMListDict(
            name="regimeMapModels", expected_type=RegimeMapModel, items=self.regimeMapModels)
        buf.write(f"{regimes_dict!r}\n")

        # Physics models
        buf.write(f"physicsModels\n{{\n")
        # Drag
        fluid1_structure_drag = OpenFOAMListDict(
            name=f"\"{self.fluid1.name}.structure\"", expected_type=DragModel, items=self.fluid1_structure.dragModels)
        fluid2_structure_drag = OpenFOAMListDict(
            name=f"\"{self.fluid2.name}.structure\"", expected_type=DragModel, items=self.fluid2_structure.dragModels)
        
        # TODO: Not sure whether ffn allows more than one fluid-fluid model (byZone?)
        if (len(self.fluid_fluid.dragModels) > 1):
            fluid_fluid_drag = OpenFOAMListDict(
                name=f"\"{self.fluid1.name}.{self.fluid2.name}\"", expected_type=DragModel, items=self.fluid_fluid.dragModels)
        else:
            fluid_fluid_drag = self.fluid_fluid.dragModels[0]
            fluid_fluid_drag.name = f"\"{self.fluid1.name}.{self.fluid2.name}\""
        
        drag_models = OpenFOAMListDict(expected_type=(OpenFOAMListDict, DragModel), name="dragModels")
        drag_models.extend([fluid1_structure_drag, fluid2_structure_drag])
        drag_models.append(fluid_fluid_drag)

        buf.write(f"{drag_models.__repr__(depth=1)}\n")

        # Two-phase drag multiplier
        buf.write(f"\ttwoPhaseDragMultiplierModel{self.fluid_fluid.twoPhaseDragMultiplierModel.__repr__(depth=1)}\n")

        # Heat transfer
        fluid1_structure_ht = OpenFOAMListDict(
            name=f"\"{self.fluid1.name}.structure\"", expected_type=HeatTransferModel, items=self.fluid1_structure.heatTransferModels)
        fluid2_structure_ht = OpenFOAMListDict(
            name=f"\"{self.fluid2.name}.structure\"", expected_type=HeatTransferModel, items=self.fluid2_structure.heatTransferModels)
        
        # TODO: Not sure whether ffn allows more than one fluid-fluid heat transfer model per fluid (byZone?)
        # There is typically at least one for fluid1 and one for fluid2
        for ht_model in self.fluid_fluid.heatTransferModels:
            ht_model.name = f"\"{ht_model.phaseName}\""
            
        fluid_fluid_ht = OpenFOAMListDict(
            name=f"\"{self.fluid1.name}.{self.fluid2.name}\"", expected_type=HeatTransferModel, items=self.fluid_fluid.heatTransferModels)
        
        ht_models = OpenFOAMListDict(expected_type=(OpenFOAMListDict, HeatTransferModel), name="heatTransferModels")
        ht_models.extend([fluid1_structure_ht, fluid2_structure_ht])
        ht_models.append(fluid_fluid_ht)

        buf.write(f"{ht_models.__repr__(depth=1)}\n")

        # Pair geometry models 
        pair_geom_models = OpenFOAMListDict(
            name="pairGeometryModels", expected_type=PairGeometryModel)
         
        if(self.fluid_fluid.pairGeometryModel is not None):
            self.fluid_fluid.pairGeometryModel.name = f"\"{self.fluid1.name}.{self.fluid2.name}\""
            pair_geom_models.append(self.fluid_fluid.pairGeometryModel)
        if(self.fluid1_structure.pairGeometryModel is not None):
            self.fluid1_structure.pairGeometryModel.name = f"\"{self.fluid1.name}.structure\""
            pair_geom_models.append(self.fluid1_structure.pairGeometryModel)    
        if(self.fluid2_structure.pairGeometryModel is not None):
            self.fluid2_structure.pairGeometryModel.name = f"\"{self.fluid2.name}.structure\""
            pair_geom_models.append(self.fluid2_structure.pairGeometryModel)
            
        buf.write(f"{pair_geom_models.__repr__(depth=1)}")        
        
        # Phase change
        buf.write(f"\tphaseChangeModel{self.fluid_fluid.phaseChangeModel.__repr__(depth=1)}\n")
        buf.write(f"}}\n\n")

        # pMin, pRefCell, pRefValue, residualKd
        buf.write(f"pMin {self.pMin};\n")
        buf.write(f"pRefCell {self.pRefCell};\n")
        buf.write(f"pRefValue {self.pRefValue};\n")
        if(self.residualKd):
            buf.write(f"residualKd {self.residualKd};\n")

    def export_to_openfoam(self):
        self.create_folders()

        self.fluid1.turbulenceProperties.region = self.region
        self.fluid1.turbulenceProperties.ext = self.fluid1.name
        self.fluid1.turbulenceProperties.export_to_openfoam()

        self.fluid1.thermophysicalProperties.region = self.region
        self.fluid1.thermophysicalProperties.ext = self.fluid1.name
        self.fluid1.thermophysicalProperties.export_to_openfoam()

        self.fluid2.turbulenceProperties.region = self.region
        self.fluid2.turbulenceProperties.ext = self.fluid2.name
        self.fluid2.turbulenceProperties.export_to_openfoam()

        self.fluid2.thermophysicalProperties.region = self.region
        self.fluid2.thermophysicalProperties.ext = self.fluid2.name
        self.fluid2.thermophysicalProperties.export_to_openfoam()

        self.export_g_to_openfoam()

        self.export_phase_properties_to_openfoam()

        self.fvSolution.extraDict['PIMPLE'] = self.pimpleOptions

        return super().export_to_openfoam()

# class CompressibleInterFoam(ThermalHydraulicsSolver):
#     """
#     OpenFOAM's solver for two compressible, non-isothermal immiscible fluids
#     using a VOF (volume of fluid) phase-fraction based interface capturing
#     approach.
#     """
#     def __init__(
#             self,
#             region="",
#             removeBaffles=False,
#             timeFolder=None,
#             mesh=None,
#             isMeshDeformation=False,
#             displacementFieldName="disp",
#             isSetFvSolutionToDefault: bool=True,
#             isSetFvSchemesToDefault: bool=True
#         ):
#         super().__init__(
#             region,
#             "compressibleInterFoam",
#             removeBaffles,
#             timeFolder,
#             mesh,
#             isMeshDeformation,
#             displacementFieldName,
#             isSetFvSolutionToDefault
#         )

#         if (isSetFvSchemesToDefault):
#             self.set_fvSchemes_default()

#         if (isSetFvSolutionToDefault):
#             self.set_fvSolution_default()


#     def set_fvSchemes_default(self):
#         self.fvSchemes = fvSchemes()

#         self.fvSchemes.ddtSchemes['default'] = 'Euler'

#         self.fvSchemes.gradSchemes['default'] = 'Gauss linear'

#         self.fvSchemes.divSchemes['div(phi_,alpha)'] = 'Gauss vanLeer'
#         self.fvSchemes.divSchemes['div(phirb,alpha)'] = 'Gauss linear'
#         for scheme in [
#             'div(rhoPhi,U)', 'div(rhoPhi,T)', 'div(rhoPhi,K)', 'div(phi,p)',
#             'div(phi,k)'
#         ]:
#             self.fvSchemes.divSchemes[scheme] = 'Gauss upwind'
#         self.fvSchemes.divSchemes['div(((rho*nuEff)*dev2(T(grad(U)))))'] = 'Gauss linear'

#         self.fvSchemes.laplacianSchemes['default'] = 'Gauss linear uncorrected'

#         self.fvSchemes.interpolationSchemes['default'] = 'linear'

#         self.fvSchemes.snGradSchemes['default'] = 'uncorrected'

#         del self.fvSchemes.fluxRequired['default']
#         self.fvSchemes.fluxRequired['alpha.water'] = None
#         self.fvSchemes.fluxRequired['p_rgh'] = None


#     def set_fvSolution_default(self):
#         self.fvSolution = fvSolution()

#         self.fvSolution.solvers['"alpha.water.*"'] = fvSolutionSolver(
#             solver="smoothSolver",
#             smoother="symGaussSeidel",
#             tolerance=1e-8,
#             relTol=0,
#             MULESCorr=True,
#             nLimiterIter=5,
#             nAlphaCorr=2,
#             nAlphaSubCycles=1,
#             cAlpha=1
#         )

#         self.fvSolution.solvers['".*(rho|rhoFinal)"'] = fvSolutionSolver(
#             solver='diagonal'
#         )

#         self.fvSolution.solvers['"pcorr.*"'] = fvSolutionSolver(
#             solver='PCG',
#             preconditioner=OpenFOAMDict({
#                 'preconditioner': 'GAMG',
#                 'tolerance': 1e-05,
#                 'relTol': 0,
#                 'smoother': 'DICGaussSeidel',
#             }),
#             tolerance=1e-5,
#             relTol=0,
#             maxIter=100
#         )

#         self.fvSolution.solvers['p_rgh'] = fvSolutionSolver(
#             solver='GAMG',
#             tolerance=1e-7,
#             relTol=0.01,
#             smoother='DIC'
#         )

#         self.fvSolution.solvers['p_rghFinal'] = fvSolutionSolver(
#             solver='PCG',
#             preconditioner=OpenFOAMDict({
#                 'preconditioner': 'GAMG',
#                 'tolerance': 1e-05,
#                 'relTol': 0,
#                 'nVcycles': 2,
#                 'smoother': 'DICGaussSeidel',
#                 'nPreSweeps': 2
#             }),
#             tolerance=1e-7,
#             relTol=0,
#             maxIter=20
#         )

#         self.fvSolution.solvers['U'] = fvSolutionSolver(
#             solver='smoothSolver',
#             smoother='symGaussSeidel',
#             tolerance=1e-6,
#             relTol=0,
#         )

#         self.fvSolution.solvers['"(T|k|B|nuTilda).*"'] = fvSolutionSolver(
#             solver='smoothSolver',
#             smoother='symGaussSeidel',
#             tolerance=1e-8,
#             relTol=0,
#         )

#         self.pimpleOptions.nOuterCorrectors = 1
#         self.pimpleOptions.nCorrectors = 3
#         self.pimpleOptions.nNonOrthogonalCorrectors = 1
#         self.pimpleOptions.momentumPredictor = False
#         self.pimpleOptions.minNOuterCorrectors = None
#         self.pimpleOptions.solveEnergy = None
#         self.pimpleOptions.solveFluidMechanics = None
#         self.pimpleOptions.momentumMode = None
#         self.pimpleOptions.correctUntilConvergence = None
#         self.pimpleOptions.porousInterfaceSharpness = None
#         self.pimpleOptions.minMagU = None


#     def export_to_openfoam(self):
#         self.create_folders()

#         self.thermophysicalProperties.region = self.region
#         self.turbulenceProperties.region = self.region
#         # self.phaseProperties.region = self.region

#         self.export_g_to_openfoam()
#         self.thermophysicalProperties.export_to_openfoam()
#         self.turbulenceProperties.export_to_openfoam()
#         # self.phaseProperties.export_to_openfoam()

#         self.fvSolution.extraDict['PIMPLE'] = self.pimpleOptions

#         # return super().export_to_openfoam()
#         return(Solver.export_to_openfoam(self))