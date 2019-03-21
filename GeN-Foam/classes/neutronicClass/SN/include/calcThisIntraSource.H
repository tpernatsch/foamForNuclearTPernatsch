
isotropicIntraSource *= 0.0;

isotropicIntraSource += flux_[energyI] * (nuSigmaEff[energyI] * chiPrompt[energyI] * (1 - BetaTot) / keff_ + sigmaFromToP[0][energyI][energyI] );

intraSources[directionI] = isotropicIntraSource;

if (xs_.legendreMoments() > 1)
{
    forAll(angularFlux_[energyI], directionJ)
    {
        volScalarField intraSourceScatteringCoeff = isotropicIntraSource*0.0; 
        for(int momentL = 1; momentL < xs_.legendreMoments(); momentL++)
        {
            intraSourceScatteringCoeff += legendreMatrices_[momentL][directionJ][directionI] * sigmaFromToP[momentL][energyI][energyI] * (2.0 * momentL + 1.0);
        }
        intraSources[directionI] += angularFlux_[energyI][directionJ] * directionWeights_[directionJ] * intraSourceScatteringCoeff ;
    }
}


