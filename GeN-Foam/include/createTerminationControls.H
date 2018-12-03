// Extra by SR

const scalar TBoil(runTime.controlDict().lookupOrDefault("TBoil", GREAT));

const bool untilSS = 
	runTime.controlDict().lookupOrDefault("untilSS", false);

int nTimeStepsSolveFluidMechanics = runTime.controlDict().lookupOrDefault("nTimeStepsSolveFluidMechanics", -1);
bool solveFluidMechanicsOnlyForFirstNTimeSteps = false;
if (nTimeStepsSolveFluidMechanics > 0)
{
	solveFluidMechanicsOnlyForFirstNTimeSteps = true;
}

int SSCounter(0);
int nTimeSteps(0);
scalar TMaxOld(0);
bool firstTimePrintTMaxt(true);