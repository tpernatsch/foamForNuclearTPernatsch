
{
	Info << "Solving ROM energy: " << energyI <<  endl;
	
	int PODBaseSize(base_fluxStar[energyI].size());

	List<scalar> sourcePOD(maxPODBaseSize,0.0);

	volScalarField Sphib("Sphib",neutroEq.SourceWithBC()); // come .source() ma volScalarField invece di scalrField. In piu' e' gia' divisa per il volume delle celle. Include BC. Dovrebbe...
	//volScalarField Sphib("Sphib",neutroEq.Source()); 

	for(int baseNoI = 0; baseNoI < PODBaseSize ; baseNoI++ )
	{
		sourcePOD[baseNoI] = Foam::POD::projection(Sphib,base_fluxStar[energyI][baseNoI]);
	}

	List<scalar> coeffPOD(maxPODBaseSize,0.0);
	solve(coeffPOD,matrixROM[energyI],sourcePOD);
	fluxStar[energyI] *= 0.0;
	for(int baseNoI = 0; baseNoI < PODBaseSize ; baseNoI++ )
	{
		fluxStar[energyI] += coeffPOD[baseNoI] * base_fluxStar[energyI][baseNoI];
	}

	fluxStar[energyI].correctBoundaryConditions();


}

