# case.py
from foamForNuclear import (
    mesh, fields, boundaryConditions as bc, offbeat_lib as offbeat, functions)
from foamForNuclear.offbeat_lib import materials
from foamForNuclear.case import OffbeatCase

# Geometry and mesh input
mm = 1e-3

fuel_length = 3000 * mm
plenum_length = 100 * mm

fuel_ri     = 0.0 * mm
fuel_ro     = 4.65 * mm

clad_ri     = 4.715 * mm
clad_ro     = 5.465 * mm

fuel_nr = 30
clad_nr = 10
slices  = 1

# default values (just for single-case runs)
lhgr = 10e3       # W/m
h_gap = 5e3       # W/m2K

time_pts_hrs = [0, 4000]

def build_case(case_name: str = "case"):
    # %% 1) Mesh
    rod_mesh = mesh.rod_1d(
        fuel_length=fuel_length,
        plenum_length=plenum_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        fuel_nr=fuel_nr,
        clad_nr=clad_nr,
        fuel_grading_nr=1,
        slices=slices,
    )

    # %% 2) Fields
    T = fields.Temperature(
        internalField=300,
        boundaryField=dict(
            fuelOuter = bc.ResistiveGap(value=300, alpha=h_gap),
            cladInner = bc.ResistiveGap(value=300, alpha=h_gap),
            cladOuter = bc.FixedValue(value=600),
        ),
    )

    # %% 3) Materials    
    uo2 = materials.Constant(name="fuel", conductivity=3)
    zircaloy = materials.Constant(name="cladding", conductivity=21)

    # %% 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=rod_mesh, caseFolder=case_name)
    
    case.thermalSolver    = offbeat.thermal_solver.SolidConduction()
    case.sliceMapper      = offbeat.slice_mapper.AutoAxialSlices()
    case.burnup           = offbeat.burnup.FromPower()

    case.heatSource = offbeat.heat_source.TimeDependentLhgr(
        timePoints=time_pts_hrs,
        lhgr=[lhgr, lhgr],
        materials=("fuel",),
    )    
    
    case.materials = [uo2, zircaloy]

    case.fields = [T]

    # %% 5) Settings + functions
    case.stressAnalysis.maxOuterIter = 100

    case.settings.endTime       = 4000
    case.settings.deltaT        = 1000
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 1
    case.settings.userTime      = "hours"


    # %% 6) FunctionObjects
    fuelCenterline = functions.Probes(
        name="fuelCenterline",
        fields=["Bu", "T"],
        probeLocations=[ [0.0, 0.0, fuel_length / 2]],
        enabled=True,
    )

    radialProfile = functions.Graph(
        name="radialProfile",
        start=[fuel_ri, 0.0, fuel_length / 2],
        end=[clad_ro, 0.0, fuel_length / 2],
        fields=["D", "Bu", "T"],
        nPoints=100,
        writeControl="timeStep",
        interpolationScheme="cellPatchConstrained"
    )

    case.functions = [fuelCenterline, radialProfile]
    
    return case, rod_mesh


if __name__ == "__main__":
    case, rod_mesh = build_case()