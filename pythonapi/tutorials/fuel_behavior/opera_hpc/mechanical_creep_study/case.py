# %% Imports

from foamForNuclear import (
    mesh, functions, fields, boundaryConditions as bc, offbeat_lib as offbeat)
from foamForNuclear.case import OffbeatCase
from foamForNuclear.common import Vector
from foamForNuclear.offbeat_lib import materials, behaviour, laws
from foamForNuclear.offbeat_lib.materials.models import ActinideDict


# %% Parameters
mm = 1e-3

fuel_length      = 110.0 * mm
plenum_length    = 0.0  * mm

clad_nr = 10
slices  = 1

# Operating conditions
time_pts_days = [0, 3000,]
p_out = 1e5
p_in = 1e7

def build_case(
    case_name: str,
    clad_ri: float,
    clad_thickness: float,
    fast_flux: float,
    T0: float
):
    # %% 1) Geometry and mesh
    clad_ro = clad_ri + clad_thickness

    rod_mesh = mesh.rod_1d(
        fuel_length=fuel_length,
        plenum_length=plenum_length,
        clad_ri=clad_ri,
        clad_ro=clad_ro,
        clad_nr=clad_nr,
        slices=slices,
        layout="clad_only"
    )

    # %% 2) Field & BCs
    T = fields.Temperature(
        internalField=T0,
        )

    D = fields.Displacement(
        internalField=[0, 0, 0],
        boundaryField=dict(
            cladInner = bc.TractionDisplacement(value=[0, 0, 0], pressure=p_in, traction=Vector(0,0,0)),
            cladOuter = bc.TractionDisplacement(value=[0, 0, 0], pressure=p_out, traction=Vector(0,0,0)),
            )
        )

    # %% 3) Materials (constant k)
    zircaloy = materials.Zircaloy(name="cladding")
    zircaloy.constitutiveLaw = laws.MisesPlasticCreep(
        yieldStress=laws.yield_stress.Constant(sigmaY=1e9),
        creep=laws.creep.Limback(relax=0.9)
    )

    # %% 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=rod_mesh, caseFolder=case_name)

    case.thermalSolver = offbeat.thermal_solver.Constant()
    case.mechanicsSolver = offbeat.mechanics_solver.SmallStrain()
    case.rheology = offbeat.rheology.PlaneStress()

    case.fastFlux = offbeat.fast_flux.TimeDependentAxialProfile(
        timePoints=time_pts_days,
        fastFlux=[fast_flux, fast_flux],
        materials=("cladding")
    )

    case.materials = [zircaloy]

    case.fields = [T, D]    


    # %%  5) Settings + functions
    case.stressAnalysis.maxOuterIter = 100

    case.settings.endTime       = 3000
    case.settings.deltaT        = 0.01
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 1
    case.settings.userTime      = "hours"

    case.settings.adjustableTimeStep         = True
    case.settings.maxDeltaT                  = 100
    case.settings.minDeltaT                  = 0.01
    case.settings.maxRelativeDeltaTIncrease  = 1e9
    case.settings.minRelativeDeltaTDecrease  = 1e9
    case.settings.maxRelativePowerIncrease   = 1e9
    case.settings.maxRelativePowerDecrease   = 1e9

    # %% 6) FunctionObjects
    cladOuter = functions.PatchProbes(
        name="cladOuter",
        fields=['D',],
        probeLocations=[[clad_ro, 0.0, fuel_length/2]],
        patches=["cladOuter",],
        enabled=True
    )

    volAverage = functions.VolFieldValue(
        name="volAverage",
        fields=["sigma"],
        operation="volAverage",
        regionName="cladding",
    )

    case.functions = [cladOuter, volAverage]

    return case, rod_mesh

if __name__ == "__main__":
    case, rod_mesh = build_case()