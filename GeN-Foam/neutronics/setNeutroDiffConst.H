Info<< "Set neutro fields\n" << endl;
{
	volScalarField DNeutroRad = DNeutro & radialOrientation;
	volScalarField DNeutroAx = DNeutro & axialOrientation;

	DNeutroRad.correctBoundaryConditions();
	DNeutroAx.correctBoundaryConditions();

	radExp = fvc::grad(DNeutroRad)& radialOrientation;
	axExp = fvc::grad(DNeutroAx)& axialOrientation;

	volScalarField logT = log(max(0.01,(TfuelNeutro/dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef))));
	volScalarField diffT = (pow(max(TfuelNeutro,dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),SMALL)),0.5)-pow(dimensionedScalar("",dimensionSet(0,0,0,1,0,0,0),TfuelRef),0.5));

	forAll(flux,energyI)
	{
		bool parametrize = true;
		forAll(doNotParametrize, groupI)
		{
			if(energyI == doNotParametrize[groupI])
			{
				parametrize = false;
				Info << "Not Parametrizing xs for group: " << doNotParametrize[groupI] << endl;
			}
		}


		forAll(entries, zoneI)
		{

			label zone = zoneI;

			const word& name = entries[zoneI].keyword();

			Info << "Mesh zone: " << name  << endl;

			label zoneId = neutroMesh.cellZones().findZoneID(name);

			if(zoneId == -1)
			{
			  	Info << "Warning: mesh zone: " << name << " does not exists" << endl;
			  	Info << "Ignore this warning if it is a completely withdrawn control rod" << endl;
			}
			else
			{


			//check if you are in a follower or CR
			bool isFollower = false;
			bool isCR = false;
			label CRzone = 0;

			forAll(CRentries,CRzoneI)
			{

				const word& CRname = CRentries[CRzoneI].keyword();
				if(name == CRname)
				{
					isCR = true ;
					CRzone = CRzoneI;
				}
				if(name == CRFollowerName[CRzoneI])
				{
					isFollower = true ;
					CRzone = CRzoneI;
				}

			}


			forAll(neutroMesh.cellZones()[zoneId], cellIlocal)
			{

				label cellIglobal = neutroMesh.cellZones()[zoneId][cellIlocal];
				zone = zoneI;

				if(CRNumber>0)
				{
					if(isFollower && neutroMesh.C()[cellIglobal].z()>CRposition[CRzone])
					{
						forAll(entries, zoneII)
						{
							if(entries[zoneII].keyword() == CRentries[CRzone].keyword())
							{
								zone = zoneII;
							}
						}
					}

					if(isCR && neutroMesh.C()[cellIglobal].z()<CRposition[CRzone])
					{
						forAll(entries, zoneII)
						{
							if(entries[zoneII].keyword() == CRFollowerName[CRzone])
							{
								zone = zoneII;
							}
						}
					}
				}



				//constant
				IV[energyI][cellIglobal] = IVList[zone][energyI];
				chiPrompt[energyI][cellIglobal] = chiPromptList[zone][energyI];
				chiDelayed[energyI][cellIglobal] = chiDelayedList[zone][energyI];
				BetaTot[cellIglobal] = BetaTotList[zone];

				//fields dependent
				D[energyI][cellIglobal] = DList[zone][energyI] + parametrize * (
									  fuelTempDList[zone][energyI]*logT[cellIglobal]*(fastNeutrons)
									+ fuelTempDList[zone][energyI]*diffT[cellIglobal]*(1-fastNeutrons)
									+ axialExpDList[zone][energyI]*(axExp[cellIglobal])
									+ radialExpDList[zone][energyI]*(radExp[cellIglobal])
									+ rhoCoolDList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
									+ TCoolDList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
									+ cladExpDList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

				nuSigmaEff[energyI][cellIglobal] = nuSigmaEffList[zone][energyI] + parametrize * (
									  fuelTempNuSigmaEffList[zone][energyI]*logT[cellIglobal]*(fastNeutrons)
									+ fuelTempNuSigmaEffList[zone][energyI]*diffT[cellIglobal]*(1-fastNeutrons)
									+ axialExpNuSigmaEffList[zone][energyI]*(axExp[cellIglobal])
									+ radialExpNuSigmaEffList[zone][energyI]*(radExp[cellIglobal])
									+ rhoCoolNuSigmaEffList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
									+ TCoolNuSigmaEffList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
									+ cladExpNuSigmaEffList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

				sigmaPow[energyI][cellIglobal] = sigmaPowList[zone][energyI] + parametrize * (
									  fuelTempSigmaPowList[zone][energyI]*logT[cellIglobal]*(fastNeutrons)
									+ fuelTempSigmaPowList[zone][energyI]*diffT[cellIglobal]*(1-fastNeutrons)
									+ axialExpSigmaPowList[zone][energyI]*(axExp[cellIglobal])
									+ radialExpSigmaPowList[zone][energyI]*(radExp[cellIglobal])
									+ rhoCoolSigmaPowList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
									+ TCoolSigmaPowList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
									+ cladExpSigmaPowList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

				sigmaDisapp[energyI][cellIglobal] = sigmaDisappList[zone][energyI] + parametrize * (
									  fuelTempSigmaDisappList[zone][energyI]*logT[cellIglobal]*(fastNeutrons)
									+ fuelTempSigmaDisappList[zone][energyI]*diffT[cellIglobal]*(1-fastNeutrons)
									+ axialExpSigmaDisappList[zone][energyI]*(axExp[cellIglobal])
									+ radialExpSigmaDisappList[zone][energyI]*(radExp[cellIglobal])
									+ rhoCoolSigmaDisappList[zone][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
									+ TCoolSigmaDisappList[zone][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
									+ cladExpSigmaDisappList[zone][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

				forAll(sigmaFromTo,energyJ)
				{
						sigmaFromTo[energyJ][energyI][cellIglobal] = sigmaFromToList[zone][energyJ][energyI] + parametrize * (
										  fuelTempSigmaFromToList[zone][energyJ][energyI]*logT[cellIglobal]*(fastNeutrons)
										+ fuelTempSigmaFromToList[zone][energyJ][energyI]*diffT[cellIglobal]*(1-fastNeutrons)
										+ axialExpSigmaFromToList[zone][energyJ][energyI]*(axExp[cellIglobal])
										+ radialExpSigmaFromToList[zone][energyJ][energyI]*(radExp[cellIglobal])
										+ rhoCoolSigmaFromToList[zone][energyJ][energyI]*(rhoCoolNeutro[cellIglobal] - rhoCoolRef)
										+ TCoolSigmaFromToList[zone][energyJ][energyI]*(TCoolNeutro[cellIglobal] - TCoolRef)
										+ cladExpSigmaFromToList[zone][energyJ][energyI]*(TcladNeutro[cellIglobal] - TcladRef) );

				}

				fuelFraction[cellIglobal] = fuelFractionList[zone];

				if(useGivenDiscFactors)
				{
					discFactor[energyI][cellIglobal] = discFactorList[zone][energyI];
				}
				else
				{
					discFactorList[zone][energyI] = discFactor[energyI][cellIglobal] ;
				}
			}//close forall cellIlocal
	
			} //close else
		}//close forall entries

		IV[energyI].correctBoundaryConditions();
		D[energyI].correctBoundaryConditions();
		nuSigmaEff[energyI].correctBoundaryConditions();
		sigmaPow[energyI].correctBoundaryConditions();
		sigmaDisapp[energyI].correctBoundaryConditions();
		forAll(sigmaFromTo,energyJ)
		{
			sigmaFromTo[energyJ][energyI].correctBoundaryConditions();
		}
		chiPrompt[energyI].correctBoundaryConditions();
		chiDelayed[energyI].correctBoundaryConditions();
		BetaTot.correctBoundaryConditions();
		fuelFraction.correctBoundaryConditions();
		discFactor[energyI].correctBoundaryConditions();	
	}
	
}
Info << "Done" << endl;

