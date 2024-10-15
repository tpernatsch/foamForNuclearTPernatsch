bool adjustTimeStep
(
    runTime.controlDict().lookupOrDefault("adjustTimeStep", false)
);

#ifndef FIRSTDELTAT
#define FIRSTDELTAT
scalar pTotOld(neutronics.power());
#endif

if (adjustTimeStep)
{
    scalar deltaT(runTime.deltaTValue());

    //- First, the thermalHydraulic class adjusts the global time step based on
    //  possible Co number limitations. The adjustTimeStep flag is read within the
    //  function itself and runTime.deltaT is set (if applicable) within it
    if (multiphysics.solveFluidMechanics() or multiphysics.solveEnergy())
    {
        deltaT = thermalHydraulics.nextDeltaT();
    }

    //- Adjusting the time step based on neutronics limitations on the first time
    //  step is useless because information relative to the "real" pTotOld is
    //  unavailable. Thus, for conservativeness, only choose the smallest time
    //  step between that specified by the thermalHydraulics (if any) and the
    //  initial one
    if (runTime.timeIndex() == runTime.startTimeIndex())
    {
        deltaT = min(deltaT, runTime.controlDict().get<scalar>("deltaT"));
    }

    //- Then, the neutronics adjusts it again, if needed. The coding of this
    //  section isn't amazing, but it is due to the fact that the thermalHydraulic
    //  class was developed independently of GeN-Foam and then shoved in it. That
    //  is why the thermalHydraulics adjusts the time step via an ad-hoc made
    //  function while the neutronics does it this way
    else if (multiphysics.solveNeutronics()) // && runTime.timeIndex() > runTime.startTimeIndex()
    {
        scalar maxDeltaT =
            runTime.controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);

        scalar maxPowerVariation =
            runTime.controlDict().lookupOrDefault<scalar>
            (
                "maxPowerVariation", 
                0.025
            );

        scalar pTot = neutronics.power();
        
        scalar powerVariation = mag((pTot - pTotOld) / (pTotOld + SMALL));
        
        pTotOld = pTot;

        scalar maxDeltaTNeutroFact = 
            mag(maxPowerVariation/(powerVariation + SMALL));

        scalar deltaTNeutroFact = 
            min(min(maxDeltaTNeutroFact, 1.0 + 0.1*maxDeltaTNeutroFact), 1.2);

        //- If solveFluidMechanics is true, the timeStep could have already
        //  been adjusted by the thermalHydraulics, in which case set the new
        //  timeStep only if the neutronics requires an even shorter time-step.
        deltaT = min
        (
            deltaTNeutroFact*deltaT,
            multiphysics.solveFluidMechanics()
                ? deltaT
                : maxDeltaT
        );
    }

    // Update deltaT once only at the end
    runTime.setDeltaT(deltaT);

    #ifdef isCommDataLayerIncluded
    // Last, after all the other setDeltaT
    if (isSolveFMI) fmu->setDeltaT();
    #endif

    Info << "deltaT = " << runTime.deltaTValue() << endl;
}
