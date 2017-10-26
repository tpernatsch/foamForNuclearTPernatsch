if (finalIter)
{
    fluidMesh.data::add("finalIteration", true);
}

if(compr)
{
    if (oCorr == 0)
    {
        #include "rhoEqn.H"
    }
}
else
{
    rho = rhoRef;
}

#include "UEqn.H"

if(finalIter)
{
    scalar residual = 0.0;
    label couplingIter = 0;
    do
    {
        Info << "Coupling iteration no:  " << couplingIter << nl << endl;
        if(solveEnergy)
        {
            #include "EEqn.H"
        }
        if(solveThermalMechanics)
        {
            #include "solveThermalMechanics.H"
            #include "displaceMeshes.H"
        }
                
        if(solveNeutronics)
        {
            #include "solveNeutronics.H"
        }
        couplingIter++;
    }while(residual > timeStepResidual && couplingIter < maxTimeStepIterations && !eigenvalueNeutronics && tighltyCoupled);

}
else
{
    scalar residual = 0.0;
    if(solveEnergy)
    {
        #include "EEqn.H"
    }
}
if(solveFluidMechanics)
{
    for (int corr=0; corr<nCorr; corr++)
    {
        #include "pEqn.H"
    }
    turb.correct();
}


if(compr)
{
    rho = thermo.rho();
}

if (finalIter)
{
    fluidMesh.data::remove("finalIteration");
    Re = porousMedium.ReynoldsNumber();
}
