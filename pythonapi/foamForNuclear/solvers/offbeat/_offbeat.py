from __future__ import annotations

from foamForNuclear.offbeat_lib import (
    materials, thermal_solver, mechanics_solver, neutronics_solver, burnup,
    fgr, gap_gas, fast_flux, heat_source, corrosion, slice_mapper,
    ElementTransportSolver, rheology, GlobalOptions, ThermoMechanicsCouplingOptions,
    StressAnalysis
)

from foamForNuclear.checkvalue import check_type
from foamForNuclear.common import *
from foamForNuclear.numerics import fvSchemes
from foamForNuclear.numerics import fvSolution, fvSolutionSolver
from foamForNuclear.mesh import Mesh
from .._solvers import Solver
from foamForNuclear.timeFolder import TimeFolder

import attrs as attr
from attrs import define, field
from foamForNuclear._attrs_tools import auto_type_validator

from io import StringIO

import foamlib

#==============================================================================*
# OFFBEAT Main Solver


@define(
    slots=True, 
    on_setattr=[attr.setters.convert, attr.setters.validate],
    field_transformer=auto_type_validator, 
    repr=False
)
class OffbeatSolver(Solver):
    """
    OFFBEAT solver parameter object.

    More info: https://gitlab.com/foam-for-nuclear/offbeat/-/tree/develop?ref_type=heads

    Parameters
    ----------
    region : str
        Name of the Offbeat region
    solver : str {"offbeat", "extendedThermoMechanics"}
        Offbeat solver name. If use in GeN-Foam, please use
        `extendedThermoMechanics`.
    thermalSolver : ThermalSolver

    mechanicsSolver : MechanicsSolver

    neutronicsSolver : NeutronicsSolver

    elementTransportSolver : ElementTransportSolver

    materials : list[materials.Material]

    rheology : Rheology

    heatSource : HeatSource

    burnup : Burnup

    fastFlux : FastFlux

    corrosion : Corrosion

    gapGasModel : GapGasModel

    fissionGasRelease : FissionGasRelease

    sliceMapper : SliceMapper

    globalOptions: GlobalOptions

    couplingOptions: ThermoMechanicsCouplingOptions

    removeBaffles : bool
        Default `False`
    timeFolder : TimeFolder
        Object representing the time folder.
    mesh : Mesh
        Mesh object
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
    couplingOptions : ThermoMechanicsCouplingOptions
    globalOptions : GlobalOptions
    thermalSolver : ThermalSolver
    mechanicsSolver : MechanicsSolver
    neutronicsSolver : NeutronicsSolver
    elementTransportSolver : ElementTransportSolver
    rheology: Rheology
    fissionGasRelease : FissionGasRelease
    gapGasModel : GapGasModel
    heatSource : HeatSource
    fastFlux : FastFlux
    corrosion : Corrosion
    burnup : Burnup
    sliceMapper : SliceMapper
    materials : OpenFOAMListDict
    stressAnalysis : StressAnalysis
    """

    solver: str = "offbeat"
    region: str = ""
    thermalSolver: thermal_solver.ThermalSolver | None = None
    mechanicsSolver: mechanics_solver.MechanicsSolver | None = None
    neutronicsSolver: neutronics_solver.NeutronicsSolver | None = None
    elementTransportSolver: ElementTransportSolver | None = None
    materials: list[materials.Material] = field(factory=list) 
    rheology: rheology.Rheology | None = None
    heatSource: heat_source.HeatSource | None = None
    burnup: burnup.Burnup | None = None
    fastFlux: fast_flux.FastFlux | None = None
    corrosion: corrosion.Corrosion | None = None
    gapGasModel: gap_gas.GapGasModel | None = None
    fissionGasRelease: fgr.FissionGasRelease | None = None
    sliceMapper: slice_mapper.SliceMapper | None = None
    globalOptions: GlobalOptions | None = field(factory=GlobalOptions)
    couplingOptions: ThermoMechanicsCouplingOptions | None = None
    removeBaffles: bool = False
    timeFolder: TimeFolder | None = None
    mesh: Mesh | None = None        
    isMeshDeformation: bool = False
    displacementFieldName: str = "disp"
    isSetFvSolutionToDefault: bool = True
    isSetFvSchemesToDefault: bool = True

    stressAnalysis: StressAnalysis = field(factory=StressAnalysis)
    
    def __attrs_post_init__(self):
        if self.isSetFvSolutionToDefault:
            self.set_fvSolution_default()
        if self.isSetFvSchemesToDefault:
            self.set_fvSchemes_default()

        super().__attrs_post_init__()

    @property
    def is_offbeat_solver(self) -> bool:
        return(self.solver == "offbeat")


    @property
    def is_extended_thermomechanics_solver(self) -> bool:
        return(self.solver == "extendedThermoMechanics")


    def add_material(self, material: materials.Material | list[materials.Material]):
        material_list = []
        if not isinstance(material, list):
            material_list.append(material)
        else:
            material_list = material

        for mat in material_list:
            check_type("material", mat, materials.Material)
            self.materials.append(mat)


    def set_fvSolution_default(self):
        self.fvSolution = fvSolution()

        # self.fvSolution.solvers["porosity"] = fvSolutionSolver(
        #     solver="PBiCG",
        #     smoother="GaussSeidel",
        #     preconditioner="DILU",
        #     tolerance=1e-10,
        #     relTol=1e-2,
        # )
        self.fvSolution.solvers["D"] = fvSolutionSolver(
            solver="PCG",
            preconditioner="FDIC",
            tolerance=1e-10,
            relTol=0.1,
        )
        self.fvSolution.solvers["T"] = copy(self.fvSolution.solvers["D"])
        self.fvSolution.solvers["neutronFlux0"] = copy(self.fvSolution.solvers["D"])

        # self.fvSolution.solvers["fuelDisp"] = fvSolutionSolver(
        #     solver="PBiCG",
        #     preconditioner="DILU",
        #     tolerance=1e-5,
        #     relTol=1e-2,
        # )
        # self.fvSolution.solvers["CRDisp"] = copy(self.fvSolution.solvers["fuelDisp"])

        self.stressAnalysis.nCorrectors = dict({"default": 1})
        self.stressAnalysis.maxOuterIter = 1000
        self.stressAnalysis.referencePairs = []
        self.stressAnalysis.D = List([1e-6, 1, 1e-6])
        self.stressAnalysis.T = 1e-6
        self.stressAnalysis.neutronFlux0 = 1e-6
        self.stressAnalysis.useRelResD = False
        self.stressAnalysis.useRelResT = False
        # self.stressAnalysis.relD = 1e-5
        # self.stressAnalysis.relT = 1e-5
        self.stressAnalysis.absErrD = 0
        self.stressAnalysis.absErrT = 0

        self.add_relaxation_on_field('D', 0.9)
        self.add_relaxation_on_field('T', 0.9)


    def set_fvSchemes_default(self):
        self.fvSchemes = fvSchemes()

        self.fvSchemes.d2dt2Schemes['default'] = 'Euler'
        self.fvSchemes.ddtSchemes['default'] = 'Euler'

        self.fvSchemes.gradSchemes['default'] = 'leastSquaresImplicitContact'

        self.fvSchemes.divSchemes['default'] = 'Gauss linear'
        for scheme in [
            'div(flux,P)', 'div(fuelDisp)', 'div(CRDisp)'
        ]:
            self.fvSchemes.divSchemes[scheme] = 'Gauss upwind'

        self.fvSchemes.laplacianSchemes['default'] = 'Gauss linear corrected'

        self.fvSchemes.interpolationSchemes['default'] = 'linear'

        self.fvSchemes.snGradSchemes['default'] = 'uncorrected'

        self.fvSchemes.fluxRequired['default'] = 'false'
        self.fvSchemes.fluxRequired['D'] = None
        self.fvSchemes.fluxRequired['T'] = None


    def get_required_fields(self) -> list[str]:
        """
        Return a list of required fields depending on the solver selected.
        They are labeled "MUST_READ"
        """
        solverFields = ["T"]
        if (self.mechanicsSolver is not None and
            self.mechanicsSolver.TYPE != "fromLatestTime" and 
            self.mechanicsSolver.TYPE != "constant"):
            solverFields += ["D"]
        if (self.is_offbeat_solver):
            solverFields += []
        if (self.is_extended_thermomechanics_solver):
            if (self.couplingOptions.correctDispForNeutro):
                solverFields += ["fuelDisp", "CRDisp"]

        return(solverFields)


    def get_default_fields(self) -> list[str]:
        """
        Return a list of default fields that can be used as input depending on
        the solver selected. They are labeled "READ_IF_PRESENT"
        """
        common = []
        solverSpecific = []

        if (self.is_offbeat_solver):
            solverSpecific = []
        elif (self.is_extended_thermomechanics_solver):
            solverSpecific = []

        return(self.get_required_fields() + solverSpecific + common)


    def import_materials_from_openfoam(
        self,
        path: str = "./"
    ) -> None:
        """
        Read materials from constant/solverDict and populate self.materials.

        Parameters
        ----------
        caseFolder : str
            Root case folder (where 'constant' lives).
        solverDictName : str
            Name of the solverDict file inside 'constant'.
        """
        # 1) open constant/solverDict with foamlib
        foam_file = foamlib.FoamFile(f"{path}/constant/solverDict")

        # foam_data should be dict-like: foam_data["materials"] → SubDict
        if "materials" not in foam_file:
            raise RuntimeError(
                f"'materials' sub-dictionary not found in {path}/constant/solverDict"
            )

        materials_dict = foam_file["materials"]

        # 2) loop over materials in the FOAM dict
        for mat_name, mat_subdict in materials_dict.items():
            # mat_name = e.g. "fuel", "cladding"
            if "type" not in mat_subdict:
                raise RuntimeError(
                    f"Material '{mat_name}' in {path}/constant/solverDict has no 'type' entry"
                )

            # 3) instantiate material; passing name helps both Python and import_from_openfoam
            mat_obj = materials.Material.import_from_openfoam(f"{path}/constant/solverDict", "materials", f"{mat_name}")
            mat_obj.name = mat_name

            # 4) add to solver
            self.add_material(mat_obj)


    def import_from_openfoam(self, path: str ='./'):
        # Import physiscs modules
        self.globalOptions = GlobalOptions.import_from_openfoam(
            f"{path}/constant/solverDict", "globalOptions")
        self.thermalSolver = thermal_solver.ThermalSolver.import_from_openfoam(
            f"{path}/constant/solverDict", "thermalSolver")
        self.mechanicsSolver = mechanics_solver.MechanicsSolver.import_from_openfoam(
            f"{path}/constant/solverDict", "mechanicsSolver")
        self.neutronicsSolver = neutronics_solver.NeutronicsSolver.import_from_openfoam(
            f"{path}/constant/solverDict", "neutronicsSolver")
        self.elementTransportSolver = ElementTransportSolver.import_from_openfoam(
            f"{path}/constant/solverDict", "elementTransport")
        self.gapGasModel = gap_gas.GapGasModel.import_from_openfoam(
            f"{path}/constant/solverDict", "gapGas")
        self.corrosion = corrosion.Corrosion.import_from_openfoam(
            f"{path}/constant/solverDict", "corrosion")
        self.heatSource = heat_source.HeatSource.import_from_openfoam(
            f"{path}/constant/solverDict", "heatSource")
        self.fast_flux = fast_flux.FastFlux.import_from_openfoam(
            f"{path}/constant/solverDict", "fastFlux")
        self.burnup = burnup.Burnup.import_from_openfoam(
            f"{path}/constant/solverDict", "burnup")
        self.fissionGasRelease = fgr.FissionGasRelease.import_from_openfoam(
            f"{path}/constant/solverDict", "fissionGasRelease")
        self.sliceMapper = slice_mapper.SliceMapper.import_from_openfoam(
            f"{path}/constant/solverDict", "sliceMapper")
        self.rheology = rheology.Rheology.import_from_openfoam(
            f"{path}/constant/solverDict", "rheology")
        
        # Import materials
        self.import_materials_from_openfoam(path)

        # Import fvSchemes and fvSolution
        self.fvSchemes.import_from_openfoam()
        self.fvSolution.import_from_openfoam()

        # Import stressAnalayis (part of fvSolution)
        self.stressAnalysis = StressAnalysis.import_from_openfoam(
            f"{path}/system/fvSolution", "stressAnalysis")


    def export_properties_to_openfoam(self):
        # pick filename
        if self.is_offbeat_solver:
            filename = "solverDict"
        elif self.is_extended_thermomechanics_solver:
            filename = "thermoMechanicalProperties"
        else:
            raise ValueError("Unknown solver type")

        # everything we might write, in the order we want it to appear
        candidates = {
            "globalOptions":     self.globalOptions,
            "thermalSolver":     self.thermalSolver,
            "mechanicsSolver":   self.mechanicsSolver,
            "neutronicsSolver":  self.neutronicsSolver,
            "elementTransport":  self.elementTransportSolver,
            "corrosion":         self.corrosion,
            "rheology":          self.rheology,
            "fissionGasRelease": self.fissionGasRelease,
            "gapGas":            self.gapGasModel,
            "heatSource":        self.heatSource,
            "fastFlux":          self.fastFlux,
            "burnup":            self.burnup,
            "sliceMapper":       self.sliceMapper,
        }

        # writer
        buf = StringIO()

        # header
        buf.write(openfoamHeader)
        buf.write(openfoamFileHeader(filename))

        if self.is_extended_thermomechanics_solver:
            buf.write(f"couplingOptions{self.couplingOptions!r}\n")
            buf.write(f"globalOptions{self.globalOptions!r}\n")
            buf.write(f"thermalSolverOptions{self.thermalSolver!r}\n")
            buf.write(f"mechanicsSolverOptions{self.mechanicsSolver!r}\n")
        elif self.is_offbeat_solver:
            # dynamic blocks (skip empties, keep order)
            for key, val in candidates.items():
                if (val is not None and not val.is_empty):
                    buf.write(f"{key}{val!r}\n")

        # materials block
        materials_dict = OpenFOAMListDict(
            name="materials", 
            expected_type=materials.Material, 
            items=self.materials
        )
        buf.write(f"{materials_dict.__repr__(depth=0)}\n")

        # footer
        buf.write(openfoamFooterLine)

        # write file
        out_path = f"constant/{self.region}/{filename}"
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(buf.getvalue())


    def export_to_openfoam(self):
        self.create_folders()

        self.export_properties_to_openfoam()

        self.fvSolution.extraDict['stressAnalysis'] = self.stressAnalysis

        return super().export_to_openfoam()
