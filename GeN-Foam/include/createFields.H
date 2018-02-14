// It is important that neutronics and solid displacemente flags are read
#include "readPhysicsToSolve.H"
#include "readSolidDisplacementFoamControls.H"

// then properties are read according to the flags
#include "readNuclearDataExt.H"
#include "readThermoMechanicalProperties.H"

// and finally the fields are created
#include "createFluidFields.H"
#include "createNeutroFields.H"
#include "createThermoMechanicalFields.H"
