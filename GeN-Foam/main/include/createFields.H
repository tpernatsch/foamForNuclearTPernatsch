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

//- Silence "unused variable" warnings, I actually do use it and have no idea
//  why the compiler does not "see" it
(void) thermalHydraulics;

//- Neutronics ---------------------------------------------------------------/

autoPtr<neutronics> pneutronics
(
    neutronics::New
    (
        neutroMesh
    )
);
neutronics& neutronics = pneutronics();

//- Silence "unused variable" warnings, I actually do use it and have no idea
//  why the compiler does not "see" it
(void) neutronics;

//- Thermo-mechanics ---------------------------------------------------------/

autoPtr<thermoMechanics> pthermoMechanics
(
    thermoMechanics::New
    (
        mechMesh
    )
);
thermoMechanics& thermoMechanics = pthermoMechanics();

//- Silence "unused variable" warnings, I actually do use it and have no idea
//  why the compiler does not "see" it
(void) thermoMechanics;
