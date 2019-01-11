// These 2 standalone flags are temporary flags:

// liquidFuel-  controls is delayed neutron precursors are transported and if
// a Schimidt number has to be given for their diffusion. It also controls if
// energy is deposited straight into the coolant. Once a fuel model is done to
// control energy deposition in a homogenous case, this flag should be split into
// a precursors control in neutronicsProperties in constant/neutroRegion and a fuel
// model inside constant/fluidRegion/fuelProperties.

// compr - activates Boussinesq behaviour on fluidMechanics solver. Upon checking
// if the Boussinesq thermophysical model is working, this flag should be eliminated
// and the fluid mechanics solver should be simplified.

const bool liquidFuel =
    runTime.controlDict().lookupOrDefault("liquidFuel", false);

const bool compr =
    runTime.controlDict().lookupOrDefault("compr", false);

// It is important that neutronics and solid displacemente flags are read
// then properties are read according to the flags
//#include "readNeutronicsProperties.H"
//#include "readNuclearDataExt.H"
//#include "readSolidDisplacementFoamControls.H"
//#include "readThermoMechanicalProperties.H"

// and finally the fields are created
#include "createFluidFields.H"
//#include "createNeutroFields.H"
//#include "createThermoMechanicalFields.H"



autoPtr<thermoMechanics> pthermoMechanics
(
    thermoMechanics::New
    (
        mechMesh
    )
);
thermoMechanics& thermoMechanics = pthermoMechanics();


autoPtr<neutronics> pneutronics
(
    neutronics::New
    (
        neutroMesh
    )
);
neutronics& neutronics = pneutronics();

SolverPerformance<vector> solvPerfU;
SolverPerformance<scalar> solvPerfp_rgh;
SolverPerformance<scalar> solvPerfE;
