# case.py
from pathlib import Path
import shutil

import foamlib
import foamForNuclear as ffn

from foamForNuclear.offbeat_lib.materials.models import ActinideDict
from foamForNuclear.case import OffbeatCase

CASE_ROOT = Path("cases/baseSweep")  # single reusable case folder
DEFAULT_T = 1000
DEFAULT_BU = 0.0

def build_case(
    T_value: float,
    Bu_value: float,
    case_name: str | None = None,
):
    """
    Build a minimal 1-cell unit cube case:
    - constant thermal solver
    - constant burnup solver
    - fixed T and Bu internal fields
    - probe at the center to extract k
    """
    # 1) Mesh: unit cube, 1 cell
    mesh = ffn.mesh.BlockMesh()    
    mesh.create_cube(
        name="fuel",
        lowX=0, lowY=0, lowZ=0,
        highX=1, highY=1, highZ=1,
        nx=1, ny=1, nz=1,
        isAddBoundaryConditions=True
    )
    
    # 2) Fields
    T = ffn.fields.Temperature(internalField=T_value)
    Bu = ffn.fields.Burnup(internalField=Bu_value)

    # 3) Materials
    fuel = ffn.offbeat_lib.materials.UO2(name="fuel",)
    fuel.isotopes = ({
        "U": ActinideDict(
            ratioOverMetal = 1.0,
            massNumbers = [235, 238],
            weightFractions = [0.045, 0.955]
    )})
    
    fuel.relocation = ffn.offbeat_lib.materials.behaviour.relocation.Relocation()
    fuel.densification = ffn.offbeat_lib.materials.behaviour.densification.Densification()
    fuel.swelling = ffn.offbeat_lib.materials.behaviour.swelling.Swelling()


    # 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=mesh, caseFolder=case_name)
    case.thermalSolver    = ffn.offbeat_lib.thermal_solver.Constant()
    case.sliceMapper      = ffn.offbeat_lib.slice_mapper.AutoAxialSlices()
    case.burnup           = ffn.offbeat_lib.burnup.Constant()

    case.materials = [fuel]
    case.fields = [T, Bu]
    case.stressAnalysis.maxOuterIter = 1000

    # 5) Settings + functions
    case.settings.application   = "offbeat"
    case.settings.endTime       = 1
    case.settings.deltaT        = 1
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 1
    case.settings.userTime      = "seconds"

    if case_name is not None:
        case.caseFolder = case_name

    # 6) FunctionObjects
    probe = ffn.functions.Probes(
        name="probeCenter",
        fields=["k"],  # extract the thermal conductivity field
        probeLocations=[[0.5, 0.5, 0.5],],
        enabled=True,
    )

    case.functions = [probe]

    return case, mesh

def clean_case_runtime(case_name: Path):
    """Remove time directories and postProcessing, keep constant/ and system/."""
    # delete numeric time dirs (e.g. 0.5, 1, 2, 10)
    for p in case_name.iterdir():
        if p.is_dir():
            try:
                float(p.name)  # numeric?
                if p.name != "0":  # keep 0/ (we will edit its fields)
                    shutil.rmtree(p)
            except ValueError:
                pass

    pp = case_name / "postProcessing"
    if pp.exists():
        shutil.rmtree(pp)

def setup_base_case():
    case, mesh = build_case(T_value=DEFAULT_T, Bu_value=DEFAULT_BU, case_name=str(CASE_ROOT))
    case.clean()                         # only once
    case.export_to_openfoam()            # only once
    ffn.run_preprocessing(case=case)             # only once
    return case, mesh

def run_case(T_value: float, Bu_value: float, case_name: str | None = None):
    """Build, clean, run, and return model."""
    case, mesh = build_case(T_value=T_value, Bu_value=Bu_value,
                              case_name=case_name)
    
    case.clean()
    case.export_to_openfoam()
    case.run()

    return case, mesh

def run_point(case: ffn.case.Case, T, Bu):
    case_name = Path(case.caseFolder) if hasattr(case, "caseFolder") else CASE_ROOT

    # clean runtime outputs, but keep mesh/dicts
    clean_case_runtime(case_name)

    # patch fields
    T_field = foamlib.FoamFieldFile(case_name / "0" / "T")
    T_field.internal_field = T
    
    Bu_field = foamlib.FoamFieldFile(case_name / "0" / "Bu")
    Bu_field.internal_field = Bu

    # run solver only
    ffn.run(case=case)

    # read probe
    probe = case.functions[0]  # probeCenter
    data, pts = probe.read_from_case(fieldName="k", startTime=0, caseFolder=str(case_name))
    return data[1][0]  # k at time=1

if __name__ == "__main__":
    # Example single test
    run_case(T_value=DEFAULT_T, Bu_value=DEFAULT_BU, case_name="test_1000K")