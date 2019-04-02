
List<PtrList<volScalarField> > sigmaFromToTempEJEI(xs_.legendreMoments());
forAll(sigmaFromToTempEJEI,momentL)
{
    sigmaFromToTempEJEI[momentL].setSize(xs_.energyGroups());
    forAll(sigmaFromToTempEJEI[momentL],energyJ)
    {
        sigmaFromToTempEJEI[momentL].set
        (
            energyJ,
            new volScalarField(xs_.sigmaFromTo(momentL,energyJ,energyI,Tfuel_,Tclad_,rhoCool_,TCool_,Disp_))
        );
    }

}

// Calc intra-group source

Info<< "Calculating intraSources" << endl;

PtrList<volScalarField> intraSources(discreteDirections_);

volScalarField isotropicIntraSource = flux_[energyI] * (nuSigmaEff[energyI] * chiPrompt[energyI] * (1 - BetaTot) / keff_ + sigmaFromToTempEJEI[0][energyI]);
//isotropicIntraSource = flux_[energyI] * (nuSigmaEff[energyI] * chiPrompt[energyI] * (1 - BetaTot) / keff_ );


forAll(angularFlux_[energyI], directionI)
{
    intraSources.set(directionI, new  volScalarField("",isotropicIntraSource));
    //intraSources[directionI] = isotropicIntraSource;
}

if (xs_.legendreMoments() > 1)
{   
    volScalarField intraSourceScatteringCoeff = nuSigmaEff[energyI]*0.0; 
    forAll(angularFlux_[energyI], directionI)
    {
        forAll(angularFlux_[energyI], directionJ)
        {
            intraSourceScatteringCoeff *= 0.0;
            for(int momentL = 1; momentL < xs_.legendreMoments(); momentL++)
            {
                scalar scatteringCoeff = legendreMatrices_[momentL][directionJ][directionI] * (2.0 * momentL + 1.0);
                intraSourceScatteringCoeff += scatteringCoeff * sigmaFromToTempEJEI[momentL][energyI] ;
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
        isotropicExtraSource += flux_[energyJ] * (nuSigmaEff[energyJ] * chiPrompt[energyI] * ((1 - BetaTot) / keff_) + sigmaFromToTempEJEI[0][energyJ] );
        
    }
}


forAll(angularFlux_[energyI], directionI)
{
    extraSources.set(directionI, new  volScalarField("", isotropicExtraSource));
}


if (xs_.legendreMoments() > 1)
{
    volScalarField extraSourceScatteringCoeff = nuSigmaEff[energyI]*0.0; 
    forAll(flux_, energyJ)
    {
        if (energyJ != energyI && xs_.sigmaFromToYesNo()[energyJ][energyI])
        {
            forAll(angularFlux_[energyI], directionI)
            {
                forAll(angularFlux_[energyJ], directionJ)
                {
                    extraSourceScatteringCoeff *= 0.0; 
                    for(int momentL = 1; momentL < xs_.legendreMoments(); momentL++)
                    {
                        scalar scatteringCoeff = legendreMatrices_[momentL][directionJ][directionI]* (2.0 * momentL + 1.0);
                        extraSourceScatteringCoeff += scatteringCoeff * sigmaFromToTempEJEI[momentL][energyJ] ;
                    }
                    extraSources[directionI] += angularFlux_[energyJ][directionJ] * extraSourceScatteringCoeff * directionWeights_[directionJ];
                }
            }
        }
    }
}


        
Info<< "Done" << endl;
