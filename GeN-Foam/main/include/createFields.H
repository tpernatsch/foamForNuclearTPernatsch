fv::options& fluidFvOptions(fv::options::New(fluidMesh));

multiphysicsControl multiphysics
(
    runTime,
    fluidMesh,
    neutroMesh,
    mechMesh
);

autoPtr<thermalHydraulicModel> pthermalHydraulics
(
    thermalHydraulicModel::New
    (
        fluidMesh, 
        multiphysics, 
        fluidFvOptions
    )
);
thermalHydraulicModel& thermalHydraulics = pthermalHydraulics();

autoPtr<neutronics> pneutronics
(
    neutronics::New
    (
        neutroMesh
    )
);
neutronics& neutronics = pneutronics();

autoPtr<thermoMechanics> pthermoMechanics
(
    thermoMechanics::New
    (
        mechMesh
    )
);
thermoMechanics& thermoMechanics = pthermoMechanics();

volScalarField& volFuelPower(thermalHydraulics.nuclearFuelPowerDensity());

#include "createNBFields.H"

/*
SolverPerformance<vector> solvPerfU;
SolverPerformance<scalar> solvPerfp_rgh;
SolverPerformance<scalar> solvPerfE;
*/