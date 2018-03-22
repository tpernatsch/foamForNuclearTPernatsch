const dictionary& neutronTransport = neutroMesh.solutionDict().subDict("neutronTransport");

const bool integralPredictor = neutronTransport.lookupOrDefault("integralPredictor", false);
const bool implicitPredictor = neutronTransport.lookupOrDefault("implicitPredictor", false);
const bool aitkenAcceleration = neutronTransport.lookupOrDefault("aitkenAcceleration", false);
const bool ROMAcceleration = neutronTransport.lookupOrDefault("ROMAcceleration", false);

scalar neutronIterationResidual =
    neutronTransport.lookupOrDefault<scalar>("neutronIterationResidual", 0.00001);

int maxNeutronIterations =
    neutronTransport.lookupOrDefault<int>("maxNeutronIterations", 10);
