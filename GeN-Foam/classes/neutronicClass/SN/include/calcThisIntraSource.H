
isotropicIntraSource *= 0.0;

isotropicIntraSource += flux_[energyI] * (nuSigmaEff[energyI] * chiPrompt[energyI] * (1 - BetaTot) / keff_ + sigmaFromToTempEIEI[0]);

intraSources[directionI] = isotropicIntraSource;

if (xs_.legendreMoments() > 1)
{
    forAll(angularFlux_[energyI], directionJ)
    {
        volScalarField intraSourceScatteringCoeff = nuSigmaEff[energyI]*0.0; 
        for(int momentL = 1; momentL < xs_.legendreMoments(); momentL++)
        {
            intraSourceScatteringCoeff += legendreMatrices_[momentL][directionJ][directionI] * sigmaFromToTempEIEI[momentL] * (2.0 * momentL + 1.0);
        }
        intraSources[directionI] += angularFlux_[energyI][directionJ] * directionWeights_[directionJ] * intraSourceScatteringCoeff ;
    }
}


