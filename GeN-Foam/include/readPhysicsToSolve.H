const bool solveFluidMechanics =
    runTime.controlDict().lookupOrDefault("solveFluidMechanics", false);

const bool solveEnergy =
    runTime.controlDict().lookupOrDefault("solveEnergy", false);

const bool solveNeutronics =
    runTime.controlDict().lookupOrDefault("solveNeutronics", false);

const bool solveThermalMechanics =
    runTime.controlDict().lookupOrDefault("solveThermalMechanics", false);
