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
            //#include "solveThermalMechanics.H"
            thermoMechanics.getFields(thermo.T(), TavFuel, TavClad, TMToFluid);
            thermoMechanics.correct(residual);
            neutronics.deformMesh(TMToNeutro,thermoMechanics.meshDisp()); 
        }
                
        if(solveNeutronics)
        {
            neutronics.getFieldsLiquidFuel(U, porosity, TavFuel, TavClad, thermo.rho()*rhok, thermo.T(), neutroToFluid, thermo, turb);
            neutronics.correct(couplingIter, residual, liquidFuel);
            volFuelPower *= 0.0;
            fluidToNeutro.mapTgtToSrc(neutronics.volFuelPower(), plusEqOp<scalar>(), volFuelPower.primitiveFieldRef());
            volFuelPower.correctBoundaryConditions();
        }
        couplingIter++;
    }while( (residual > timeStepResidual) && (couplingIter < maxTimeStepIterations) && !neutronics.eigenvalue() && tightlyCoupled);

}
else
{
    if(solveFluidMechanics)//first pressure corrector before energy to avoid that non-continuity of mass generate spurious heat sources 
    {                       //(can be a probloem with large time steps)
        int corr = 0;
        #include "pEqn.H"
    }
            
    scalar residual = 0.0;
    if(solveEnergy)
    {
        #include "EEqn.H"
    }
}
if(solveFluidMechanics)
{
    for (int corr=1; corr<nCorr; corr++)//start from 1 because first corrector done above
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
