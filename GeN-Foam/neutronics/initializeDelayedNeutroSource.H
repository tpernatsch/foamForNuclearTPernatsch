
delayedNeutroSource *= 0.0;
forAll(prec,precI)
{
	delayedNeutroSource += lambda[precI]*prec[precI];
}


