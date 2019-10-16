//- As of now, residuals passed to thermalHydraulics (both
//  correctFluidMechanics and correctEnergy) serve no purpose

scalar residual = 0.0;
if (multiphysics.solveFluidMechanics())
{
    thermalHydraulics.correctRegimes();
    thermalHydraulics.correctFluidMechanics(residual);
}

if(!multiphysics.finalIter())
{
    if(multiphysics.solveEnergy())
    {
        thermalHydraulics.correctEnergy(residual);
    }
}
else
{
    scalar residual = 0.0;
    label couplingIter = 0;
    do
    {
        Info << "Coupling iteration no:  " << couplingIter << nl << endl;
        if(multiphysics.solveEnergy())
        {
            thermalHydraulics.correctEnergy(residual);
        }
        #include "transferFieldsNB.H"
        if(multiphysics.solveThermoMechanics())
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
        if(multiphysics.solveNeutronics())
        {
            neutronics.getFields
            (
                TavFuelNB, 
                TavCladNB, 
                rhoNB, 
                TNB, 
                neutroToFluid
            );
            neutronics.correct(couplingIter, residual, false);
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
            (residual > multiphysics.timeStepResidual())
        &&  (couplingIter < multiphysics.maxTimeStepIterations()) 
        &&  (!neutronics.eigenvalue() && multiphysics.tightlyCoupled())
        );
}
