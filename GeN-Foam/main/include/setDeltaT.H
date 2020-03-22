//- Declare and initialize pTotOld only the first time setDeltaT.H is included
#ifndef PTOTOLD
#define PTOTOLD
scalar pTotOld(0.0);
#endif

//- First, the thermalHydraulic class adjusts the global time step based on
//  possible Co number limitations. The adjustTimeStep flag is read within the
//  fuction itself (and then read against outside, for the neutronics... read
//  comment below)
if (multiphysics.solveFluidMechanics())
{
    thermalHydraulics.adjustTimeStep();
}

//- Then, the neutronics adjusts it again, if needed. The coding of this
//  section isn't amazing, but it is due to the fact that the thermalHydraulic
//  class was developed independently of GeN-Foam and then shoved in it. That
//  is why the thermalHydraulics adjusts the time step via an ad-hoc made
//  fuction while the neutronics does it this way
if (runTime.controlDict().lookupOrDefault("adjustTimeStep", false))
{
    scalar maxDeltaT =
        runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

    scalar maxPowerVariation =
        runTime.controlDict().lookupOrDefault<scalar>
        (
            "maxPowerVariation", 
            0.025
        );

    if (multiphysics.solveNeutronics())
    {
        scalar pTot = neutronics.power();
        
        scalar powerVariation = mag((pTot - pTotOld) / (pTotOld + SMALL));
        
        pTotOld = pTot;

        scalar maxDeltaTNeutroFact = 
            mag(maxPowerVariation/(powerVariation + SMALL));

        scalar deltaTNeutroFact = 
            min(min(maxDeltaTNeutroFact, 1.0 + 0.1*maxDeltaTNeutroFact), 1.2);

        //- If solveFluidMechanics is true, the timeStep could have already
        //  been adjusted by the thermalHyraulics, in which case set the new
        //  timeStep only if the neutronics requires an even shorter time-step.
        runTime.setDeltaT
        (
            min
            (
                deltaTNeutroFact*runTime.deltaTValue(),
                (multiphysics.solveFluidMechanics()) ?
                runTime.deltaTValue() : maxDeltaT
            )
        );
    }

    Info << "deltaT = " << runTime.deltaTValue() << endl;
}
    