neutroResidual = 0.0;
/*
if(aitkenAcceleration)
{
    #include "aitkenAcceleration.H"
}   
*/

forAll(flux_,energyI)
{
    #include "calcScatteringSourceSN.H"
    scalar innerResidual = 0.0;
    int innerIter = 0;

    do
    {
        innerResidual = 0.0;
        forAll(angularFlux_[energyI], directionI)
        {
            #include "calcThisIntraSource.H"
            scalar localResidual =  solve
            (
                neutroEqList[energyI][directionI]
                ==
                extraSources[directionI]
                + intraSources[directionI]
                + delayedNeutroSource_*chiDelayed[energyI]
            ).initialResidual();
        innerResidual = max(innerResidual, localResidual);
        neutroResidual = max(neutroResidual, localResidual);
        if(eigenvalueNeutronics_) angularFlux_[energyI][directionI] = max(angularFlux_[energyI][directionI],(angularFlux_[energyI][directionI]*0.0));
        #include "updateThisScalarFlux.H"
        }
        innerIter++;
    }
    while((innerResidual > (neutroResidual*innToOutTolerance))&&(innerIter < maxInnerIter));

}




