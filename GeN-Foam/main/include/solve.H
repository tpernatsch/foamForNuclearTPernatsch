//- As of now, residuals passed to thermalHydraulics (both
//  correctFluidMechanics and correctEnergy) serve no purpose
scalar residual = 0.0;

bool solveFM(multiphysics.solveFluidMechanics());
bool solveE(multiphysics.solveEnergy());
bool solveN(multiphysics.solveNeutronics());
bool solveTM(multiphysics.solveThermoMechanics());

//- Correct flow regime map
thermalHydraulics.correctRegimes(solveFM, solveE);

if (solveFM) thermalHydraulics.correctFluidMechanics(residual);

if (!multiphysics.finalIter())
{
    if (solveE) thermalHydraulics.correctEnergy(residual);
}
else if (solveE or solveN or solveTM)
{
    scalar residual = 0.0;
    label couplingIter = 0;
    do
    {
        Info << "Coupling: iteration " << couplingIter << nl << endl;
        if (solveE)
        {
            thermalHydraulics.correctEnergy(residual);
        }
        #include "transferFieldsNB.H"
        if (solveTM)
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
        if (solveN)
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
            (residual > multiphysics.timeStepResidual()) //- Currently useless
        &&  (couplingIter < multiphysics.maxTimeStepIterations()) 
        &&  (!neutronics.eigenvalue() && multiphysics.tightlyCoupled())
        );
}
