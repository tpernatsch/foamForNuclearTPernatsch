
neutroSource *= 0.0;
forAll(flux,energyI)
{
	neutroSource += nuSigmaEff[energyI]*flux[energyI];
}

