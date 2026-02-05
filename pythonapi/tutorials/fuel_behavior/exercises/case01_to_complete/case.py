# %% Imports

from foamForNuclear import (
    mesh, functions, fields, boundaryConditions as bc, offbeat_lib as offbeat)
from foamForNuclear.case import OffbeatCase
from foamForNuclear.offbeat_lib import materials



# %% Parameters (edit here)
mm = 1e-3

fuel_length      = 110.0 * mm
plenum_length    = 10.0  * mm

# TODO: complete
fuel_ri          = 
fuel_ro          = 

gap              = 
clad_thickness   = 

clad_ri          = fuel_ro + gap
clad_ro          = clad_ri + clad_thickness

fuel_nr = 
clad_nr = 
slices  = 

# TODO: complete
T0    =       # K
Tcold =       # K
h_gap =       # W/(m^2 K)

# TODO: complete
# Power profile (days, W/m)
time_pts_days = []
lhgr_W_per_m  = []

def build_case():
    # %% 1) Geometry and mesh
    rod_mesh = mesh.rod_1d(
        fuel_length=fuel_length,
        plenum_length=plenum_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        fuel_nr=fuel_nr,
        clad_nr=clad_nr,
        slices=slices,
    )

    # %% 2) Field & BCs
    T = fields.Temperature(
        internalField=T0,
        boundaryField=dict(
            fuelOuter = bc.ResistiveGap(value=T0, alpha=h_gap),
            cladInner = bc.ResistiveGap(value=T0, alpha=h_gap),
            cladOuter = bc.FixedValue(value=Tcold),
            )
        )

    # %% 3) Materials (constant k)
    # TODO: complete
    uo2 = materials.Constant(name="fuel", conductivity=)
    zircaloy = materials.Constant()

    # %% 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=rod_mesh, caseFolder="case")

    case.thermalSolver = offbeat.thermal_solver.SolidConduction()
    case.heatSource = offbeat.heat_source.TimeDependentLhgr(
        timePoints=time_pts_days,
        lhgr=lhgr_W_per_m,
        materials=("fuel",)
    )

    case.materials = [uo2, zircaloy]
    case.fields = [T]    

    # %%  5) Settings + functions
    case.stressAnalysis.maxOuterIter = 100

    case.settings.endTime       = 365
    case.settings.deltaT        = 0.01
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 1
    case.settings.userTime      = "days"

    case.settings.adjustableTimeStep         = True
    case.settings.maxDeltaT                  = 7
    case.settings.minDeltaT                  = 0.01
    case.settings.maxRelativeDeltaTIncrease  = 1e9
    case.settings.minRelativeDeltaTDecrease  = 1e9
    case.settings.maxRelativePowerIncrease   = 1e9
    case.settings.maxRelativePowerDecrease   = 1e9

    # %% 6) FunctionObjects
    fuelCenterline = functions.Probes(
        name="fuelCenterline",
        fields=['T'],
        probeLocations=[[fuel_ri, 0.0, fuel_length/2]],
        enabled=True
    )

    radialProfile = functions.Graph(
        name="radialProfile",
        start=[fuel_ri, 0.0, fuel_length/2],
        end=[fuel_ro, 0.0,  fuel_length/2],
        nPoints=30,
        fields=['T'],
        interpolationScheme="cell"
    )

    case.functions = [fuelCenterline, radialProfile]

    return case, rod_mesh

if __name__ == "__main__":
    case, rod_mesh = build_case()