

PtrList<volScalarField> intraSources(discreteDirections_);

volScalarField isotropicIntraSource = flux_[0] * (nuSigmaEff[0] * chiPrompt[0] * (1 - BetaTot) / keff_ );


forAll(angularFlux_,energyI)
{
	forAll(angularFlux_[energyI], directionI)
	{
	    intraSources.set(directionI, new  volScalarField(isotropicIntraSource));
	}
}

PtrList<volScalarField> extraSources(discreteDirections_);

volScalarField isotropicExtraSource = isotropicIntraSource * 0.0;

forAll(angularFlux_,energyI)
{
	forAll(angularFlux_[energyI], directionI)
	{
	    extraSources.set(directionI, new  volScalarField(isotropicExtraSource));
	}
}


