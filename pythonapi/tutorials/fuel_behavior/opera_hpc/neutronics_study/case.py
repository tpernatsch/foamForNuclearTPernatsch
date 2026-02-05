# case.py
from foamForNuclear import (
    mesh, fields, boundaryConditions as bc, offbeat_lib as offbeat, functions)
from foamForNuclear.offbeat_lib import materials
from foamForNuclear.case import OffbeatCase
from foamForNuclear.offbeat_lib.materials.models import ActinideDict

# Geometry and mesh input
mm = 1e-3
fuel_length = 1.0 * mm
fuel_ri = 0.0 * mm
fuel_ro = 4.65 * mm

fuel_nr = 1000
slices = 1

# Rod history
time_pts_days = [0, 1096]
lhgr_W_per_m = [60e3, 60e3]


def build_case(
        case_name: str,
        enrichment: float
    ):
    # %% 1) Mesh
    rod_mesh = mesh.rod_1d(
        fuel_length=fuel_length,
        fuel_ri=fuel_ri,
        fuel_ro=fuel_ro,
        fuel_nr=fuel_nr,
        fuel_grading_nr=0.1,
        slices=slices,
        layout="fuel_only",
    )

    # %% 2) Fields
    neutron_flux = fields.NeutronFlux0(
        internalField=0.0,
        boundaryField={
            "fuelOuter": bc.FixedValue(value=1),
        },
    )

    # %% 3) Materials
    uo2 = materials.UO2(name="fuel")
    uo2.isotopes = ({
        "U": ActinideDict(
            ratioOverMetal = 1.0,
            massNumbers = [235, 238],
            weightFractions = [enrichment, 1 - enrichment]
    )})

    # %% 4) OffbeatCase and physics setup
    case = OffbeatCase(mesh=rod_mesh, caseFolder=case_name)

    case.neutronicsSolver = offbeat.neutronics_solver.Diffusion()
    case.burnup = offbeat.burnup.Lassmann(convergencePrecision=1e-6)

    case.heatSource = offbeat.heat_source.TimeDependentLhgr(
        timePoints=[0, 1096],
        lhgr=[60e3, 60e3],
        materials=("fuel",),
    )

    case.materials = [uo2]

    case.fields = [neutron_flux]

    # %% 5) Settings + functions
    case.stressAnalysis.maxOuterIter = 1

    case.settings.endTime       = 600
    case.settings.deltaT        = 1
    case.settings.writeControl  = "timeStep"
    case.settings.writeInterval = 100
    case.settings.userTime      = "days"

    # %% 6) FunctionObjects
    fuelCenterline = functions.Probes(
        name="fuelCenterline",
        fields=["Bu", "T"],
        probeLocations=[[0.0, 0.0, fuel_length / 2]],
        enabled=True,
    )

    radialProfile = functions.Graph(
        name="radialProfile",
        start=[fuel_ri, 0.0, fuel_length / 2],
        end=[fuel_ro, 0.0, fuel_length / 2],
        fields=["D", "Bu", "T", "N_U", "N_U235", "N_U238",
                "N_Pu238", "N_Pu239", "N_Pu240", "N_Pu241", "N_Pu242"],
        nPoints=fuel_nr,
        writeControl="timeStep",
    )

    volAverage = functions.VolFieldValue(
        name="volAverage",
        fields=["D", "Bu", "T", "N_U", "N_U235", "N_U238",
                "N_Pu238", "N_Pu239", "N_Pu240", "N_Pu241", "N_Pu242"],
        operation="volAverage",
        regionName="fuel",
    )

    sumFields = functions.VolFieldValue(
        name="sum",
        fields=["N_U", "N_U235", "N_U238",
                "N_Pu238", "N_Pu239", "N_Pu240", "N_Pu241", "N_Pu242"],
        operation="volIntegrate",
        regionName="fuel",
    )

    case.functions = [fuelCenterline, radialProfile, volAverage, sumFields]
    
    return case, rod_mesh


if __name__ == "__main__":
    case, rod_mesh = build_case()