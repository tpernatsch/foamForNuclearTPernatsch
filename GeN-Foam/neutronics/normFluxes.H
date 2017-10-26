
forAll(fluxStar, energyI)
{
    fluxStar[energyI] *= pTarget/pTot;
    fluxStar2[energyI] *= pTarget/pTot;
    flux[energyI]=(fluxStar[energyI]-fluxStar2[energyI]*2.0*SP3Neutronics)/discFactor[energyI];
}

