
forAll(fluxStar, energyI)
{
    fluxStar[energyI] *= pTarget_/pTot;
    flux[energyI]=fluxStar[energyI]/discFactor[energyI];
}

