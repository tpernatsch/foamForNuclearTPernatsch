//- FvOptions
fv::options& fluidFvOptions(fv::options::New(fluidMesh));

//- Multiphysics solution control
multiphysicsControl multiphysics
(
    runTime,
    fluidMesh,
    neutroMesh,
    mechMesh
);

//- Thermal-hydraulics -------------------------------------------------------/

autoPtr<thermalHydraulicModel> pthermalHydraulics
(
    thermalHydraulicModel::New
    (
        runTime,
        fluidMesh, 
        multiphysics, 
        fluidFvOptions
    )
);
thermalHydraulicModel& thermalHydraulics = pthermalHydraulics();

volScalarField dummyPowerDensity
(
    IOobject
    (
        "dummyPowerDensity",
        fluidMesh.time().timeName(),
        fluidMesh.time()
    ),
    fluidMesh,
    dimensionedScalar("", dimPower/dimVol, 0.0)
);
volScalarField& volFuelPower(dummyPowerDensity);

//- Neutronics ---------------------------------------------------------------/

autoPtr<neutronics> pneutronics
(
    neutronics::New
    (
        neutroMesh
    )
);
neutronics& neutronics = pneutronics();

//- Thermo-mechanics ---------------------------------------------------------/

autoPtr<thermoMechanics> pthermoMechanics
(
    thermoMechanics::New
    (
        mechMesh
    )
);
thermoMechanics& thermoMechanics = pthermoMechanics();

#include "createNBFields.H"
