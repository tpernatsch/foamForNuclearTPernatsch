// We do not have a top-level mesh. Construct the fvSolution for
// the runTime instead.
fvSolution solutionDict(runTime);

const dictionary& pimple = solutionDict.subDict("PIMPLE");

const int nOuterCorr =
    pimple.lookupOrDefault<int>("nOuterCorrectors", 1);

const bool tighltyCoupled =
    pimple.lookupOrDefault("tighltyCoupled", true);

scalar timeStepResidual =
    pimple.lookupOrDefault<scalar>("timeStepResidual", 0.001);

int maxTimeStepIterations =
    pimple.lookupOrDefault<int>("maxTimeStepIterations", 5);
