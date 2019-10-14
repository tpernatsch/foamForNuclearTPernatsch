if(compr)
{
    if (multiphysics.firstIter())
    {
        #include "rhoEqn.H"
    }
}
else
{
    rho = rhoRef;
}

#include "UEqn.H"

if(!multiphysics.finalIter())
{
    if (solveFluidMechanics)
    {
        while(multiphysics.correct())
        {
            #include "pEqn.H"
        } 
    }
    scalar residual = 0.0;
    if(solveEnergy)
    {
        #include "EEqn.H"
    }
    Re = porousMedium.ReynoldsNumber();
}
else
{
    if (solveFluidMechanics)
    {
        while(multiphysics.correct())
        {
            #include "pEqn.H"
        }
    }
    scalar residual = 0.0;
    label couplingIter = 0;
    do
    {
        Info << "Coupling iteration no:  " << couplingIter << nl << endl;
        if(solveEnergy)
        {
            #include "EEqn.H"
        }
        #include "transferFieldsNB.H"
        if(solveThermalMechanics)
        {
            thermoMechanics.getFields
            (
                TstructuresNB, 
                TavFuelNB, 
                TavCladNB, 
                TMToFluid
            );
            thermoMechanics.correct(residual);
            neutronics.deformMesh(TMToNeutro,thermoMechanics.meshDisp()); 
        }       
        if(solveNeutronics)
        {
            neutronics.getFieldsLiquidFuel
            (
                UNB, 
                porosityNB, 
                TavFuelNB, 
                TavCladNB, 
                rhoNB*rhokNB, 
                TNB, 
                muNB, 
                alphatNB,  
                neutroToFluid
            );
            neutronics.correct(couplingIter, residual, liquidFuel);
            volFuelPower *= 0.0;
            fluidToNeutro.mapTgtToSrc
            (
                neutronics.volFuelPower(), 
                plusEqOp<scalar>(), 
                volFuelPower.primitiveFieldRef()
            );
            volFuelPower.correctBoundaryConditions();
        }
        couplingIter++;
    }   
    while
        ( 
            (residual > timeStepResidual)
        &&  (couplingIter < maxTimeStepIterations) 
        &&  !neutronics.eigenvalue() && tightlyCoupled
        );
    Re = porousMedium.ReynoldsNumber();
}

if(compr)
{
    rho = thermo.rho();
}
