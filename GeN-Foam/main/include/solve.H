//- As of now, this serves no purpose
scalar residual = 0.0;

bool solveFM(multiphysics.solveFluidMechanics());
bool solveE(multiphysics.solveEnergy());
bool solveN(multiphysics.solveNeutronics());
bool solveTM(multiphysics.solveThermoMechanics());

//- Correct flow regime map
thermalHydraulics.correctRegimes(solveFM, solveE);

//- Solve fluid-mechanics
if (solveFM) thermalHydraulics.correctFluidMechanics(residual);

//- Solve energy
if (!multiphysics.finalIter())
{
    if (solveE) thermalHydraulics.correctEnergy(residual);
}

//- Solve energy-neutronics-thermomechanics coupling on last iteration
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
        if (solveTM or solveN)
        {
            #include "correctCouplingFields.H"
        }
        if (solveTM)
        {
            thermoMechanics.getFields
            (
                bafflelessTStruct, 
                bafflelessTFuelAv, 
                bafflelessTCladAv, 
                TMToFluid
            );
            thermoMechanics.correct(residual);
            neutronics.deformMesh(TMToNeutro,thermoMechanics.meshDisp()); 
        }       
        if (solveN)
        {
            neutronics.getFieldsLiquidFuel
            (
                bafflelessU,
                bafflelessAlpha,
                bafflelessTFuelAv, 
                bafflelessTCladAv, 
                bafflelessRho, 
                bafflelessTCool,
                bafflelessNullField,
                bafflelessNullField,
                neutroToFluid
            );
            neutronics.correct
            (
                couplingIter, 
                residual, 
                multiphysics.liquidFuel()
            );
            (*volFuelPower) *= 0.0;
            fluidToNeutro.mapTgtToSrc
            (
                neutronics.volFuelPower(), 
                plusEqOp<scalar>(), 
                volFuelPower->primitiveFieldRef()
            );
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
