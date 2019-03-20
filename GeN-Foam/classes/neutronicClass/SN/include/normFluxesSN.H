
forAll(flux_, energyI)
{
	forAll(angularFlux_[energyI], directionI)
	{
		angularFlux_[energyI][directionI] *= pTarget_/pTot;
	}
}
