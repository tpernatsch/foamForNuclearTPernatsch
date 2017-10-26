const bool tighltyCoupled =
    runTime.controlDict().lookupOrDefault("tighltyCoupled", true);

const bool solveFluidMechanics =
    runTime.controlDict().lookupOrDefault("solveFluidMechanics", false);

const bool solveEnergy =
    runTime.controlDict().lookupOrDefault("solveEnergy", false);

const bool solveNeutronics =
    runTime.controlDict().lookupOrDefault("solveNeutronics", false);

const bool solveThermalMechanics =
    runTime.controlDict().lookupOrDefault("solveThermalMechanics", false);

const bool eigenvalueNeutronics =
    runTime.controlDict().lookupOrDefault("eigenvalueNeutronics", false);

const bool SP3Neutronics =
    runTime.controlDict().lookupOrDefault("SP3Neutronics", false);

const bool fastNeutrons =
    runTime.controlDict().lookupOrDefault("fastNeutrons", false);

const bool liquidFuel =
    runTime.controlDict().lookupOrDefault("liquidFuel", false);

const bool compr =
    runTime.controlDict().lookupOrDefault("compr", false);

const bool adjustDiscFactors =
    runTime.controlDict().lookupOrDefault("adjustDiscFactors", false);

const bool useGivenDiscFactors =
    runTime.controlDict().lookupOrDefault("useGivenDiscFactors", false);

List<int> groupsWoDF(runTime.controlDict().lookupOrDefault<List<int> >("groupsWoDF", List<int>()));

const bool integralPredictor =
    runTime.controlDict().lookupOrDefault("integralPredictor", false);

const bool implicitPredictor =
    runTime.controlDict().lookupOrDefault("implicitPredictor", false);

const bool aitkenAcceleration =
    runTime.controlDict().lookupOrDefault("aitkenAcceleration", false);

const bool ROMAcceleration = runTime.controlDict().lookupOrDefault("ROMAcceleration", false);

List<int> doNotParametrize(runTime.controlDict().lookupOrDefault<List<int> >("doNotParametrize", List<int>()));

