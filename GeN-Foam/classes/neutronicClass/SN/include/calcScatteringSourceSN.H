
// Calc intra-group source

Info<< "Calculating intraSources" << endl;

PtrList<volScalarField> intraSources(discreteDirections_);

volScalarField isotropicIntraSource = flux_[energyI] * (nuSigmaEff[energyI] * chiPrompt[energyI] * (1 - BetaTot) / keff_ + sigmaFromToP[0][energyI][energyI] );

forAll(angularFlux_[energyI], directionI)
{
    intraSources.set(directionI, new  volScalarField(isotropicIntraSource));
}

if (xs_.legendreMoments() > 1)
{
    forAll(angularFlux_[energyI], directionI)
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
}
        
Info<< "Done" << endl;


// Calc extra-group source

Info<< "Calculating extraSources" << endl;

PtrList<volScalarField> extraSources(discreteDirections_);

volScalarField isotropicExtraSource = isotropicIntraSource * 0.0;

forAll(flux_, energyJ)
{
    if (energyJ != energyI)
    {
        isotropicExtraSource += flux_[energyJ] * (nuSigmaEff[energyJ] * chiPrompt[energyI] * (1 - BetaTot) / keff_ + sigmaFromToP[0][energyJ][energyI] );
        
    }
}


forAll(angularFlux_[energyI], directionI)
{
    extraSources.set(directionI, new  volScalarField(isotropicExtraSource));
}

if (xs_.legendreMoments() > 1)
{
    forAll(flux_, energyJ)
    {
        if (energyJ != energyI)//&&(sigmaFromToP[0][energyJ][energyI] > 0.0))
        {
            forAll(angularFlux_[energyI], directionI)
            {
                forAll(angularFlux_[energyJ], directionJ)
                {
                    volScalarField extraSourceScatteringCoeff = isotropicExtraSource*0.0; 
                    for(int momentL = 1; momentL < xs_.legendreMoments(); momentL++)
                    {
                        extraSourceScatteringCoeff += legendreMatrices_[momentL][directionJ][directionI] * sigmaFromToP[momentL][energyJ][energyI] * (2.0 * momentL + 1.0);
                    }
                    extraSources[directionI] += angularFlux_[energyJ][directionJ] * directionWeights_[directionJ] * extraSourceScatteringCoeff ;
                }
            }
        }
    }
}
        
Info<< "Done" << endl;
