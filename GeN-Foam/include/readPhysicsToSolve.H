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

const bool tighltyCoupled =
    runTime.controlDict().lookupOrDefault("tighltyCoupled", true);

const bool liquidFuel =
    runTime.controlDict().lookupOrDefault("liquidFuel", false);

const bool adjustDiscFactors =
    runTime.controlDict().lookupOrDefault("adjustDiscFactors", false);

const bool useGivenDiscFactors =
    runTime.controlDict().lookupOrDefault("useGivenDiscFactors", false);
