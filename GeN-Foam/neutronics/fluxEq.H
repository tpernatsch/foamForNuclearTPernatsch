{

	if(aitkenAcceleration)
	{
		#include "aitkenAcceleration.H"
	}	

	#include "solveROMConditions.H"

	forAll(flux,energyI)
	{
		Dalbedo = D[energyI];
		fluxStarAlbedo = fluxStar2[energyI]*(SP3Neutronics);
		#include "calcScatteringSource.H"
		neutroSource -= nuSigmaEff[energyI]*flux[energyI]*(1-eigenvalueNeutronics);

		fvMatrixExt<scalar> neutroEq = neutroEqList[energyI];


		if(!SP3Neutronics)
		{
			neutroEq = neutroEq
			  - neutroSource*(1.0-BetaTot)*chiPrompt[energyI]/keff
			  - delayedNeutroSource*chiDelayed[energyI]
			  - scatteringSourceExtra ;
		}
		else
		{
			neutroEq = neutroEq
			  - neutroSource*(1.0-BetaTot)*chiPrompt[energyI]/keff
			  - delayedNeutroSource*chiDelayed[energyI]
			  - scatteringSourceExtra
			  - (sigmaDisapp[energyI]*(fluxStar2[energyI]/discFactor[energyI])*2.0) 
			  - fvc::ddt((1-eigenvalueNeutronics)*IV[energyI]/discFactor[energyI],fluxStar2[energyI])*2.0 ;
		}

		scalar initRes = neutroEq.initialResidual();
		neutroResidual = max(neutroResidual,initRes);

		Info<< "neutro Residual = " << initRes << nl << endl;

		if(ROMAcceleration)
		{

			if( solveROM )
			{
				#include "ROMSolution.H" 
			}
			else
			{
				Info<< "ROM verification step"  << nl << endl;
				SolverPerformance<scalar> neutroSolverPerformance(neutroEq.solve());
				#include "adjustROMBase.H"
			}
	
		}
		else
		{
			Info<< "No ROM acc" << nl << endl;
			neutroEq.solve();
		}

	
		if(SP3Neutronics)
		{	
			fluxStarAlbedo = fluxStar[energyI];

			/*
			Dalbedo = D[energyI]/discFactor[energyI]*9.0/7.0; //like in PARCS
			fvScalarMatrix neutroEq2 //like in PARCS
			(
			  fvm::ddt((1-eigenvalueNeutronics)*IV[energyI]/discFactor[energyI]*3.0,fluxStar2[energyI])
			  - fvm::laplacian((D[energyI]/discFactor[energyI]*9.0/7.0),fluxStar2[energyI])
			  + fvm::Sp(((1/D[energyI]/3.0)*5.0/3.0+(sigmaDisapp[energyI]*4.0/3.0)),fluxStar2[energyI]) 
			  - (sigmaDisapp[energyI]*(fluxStar[energyI]/discFactor[energyI]))*2.0/3.0
			  + neutroSource*(1.0-BetaTot)*chiPrompt[energyI]/keff*2.0/3.0
			  + delayedNeutroSource*chiDelayed[energyI]*2.0/3.0
			  + scatteringSourceExtra*2.0/3.0
			  - fvc::ddt((1-eigenvalueNeutronics)*IV[energyI]/discFactor[energyI],fluxStar[energyI])*2.0/3.0
			 );
			*/
			Dalbedo = 1.0/(sigmaDisapp[energyI]+sigmaFromTo[energyI][energyI])*3.0/7.0;
			fvScalarMatrix neutroEq2
			(
			  fvm::ddt((1-eigenvalueNeutronics)*IV[energyI]/discFactor[energyI]*3.0,fluxStar2[energyI])
			  - fvm::laplacian((discFactor[energyI]/(sigmaDisapp[energyI]+sigmaFromTo[energyI][energyI])*3.0/7.0),fluxStar2[energyI])
			  + fvm::Sp(((sigmaDisapp[energyI]+sigmaFromTo[energyI][energyI])*5.0/3.0+(sigmaDisapp[energyI]*4.0/3.0)),fluxStar2[energyI]) 
			  - (sigmaDisapp[energyI]*(fluxStar[energyI]/discFactor[energyI]))*2.0/3.0
			  + neutroSource*(1.0-BetaTot)*chiPrompt[energyI]/keff*2.0/3.0
			  + delayedNeutroSource*chiDelayed[energyI]*2.0/3.0
			  + scatteringSourceExtra*2.0/3.0
			  - fvc::ddt((1-eigenvalueNeutronics)*IV[energyI]/discFactor[energyI],fluxStar[energyI])*2.0/3.0
			 );

			neutroEq2.solve();
		}


		flux[energyI]=(fluxStar[energyI]-fluxStar2[energyI]*2.0*SP3Neutronics)/discFactor[energyI];
		flux[energyI].correctBoundaryConditions();
		neutroSource += nuSigmaEff[energyI]*flux[energyI]*(1-eigenvalueNeutronics);

		volFuelPowerNeutro += flux[energyI] * sigmaPow[energyI];


	}
	neutroResidual0 = neutroResidual1;
	neutroResidual1 = neutroResidual2;
	neutroResidual2 = neutroResidual;

}

