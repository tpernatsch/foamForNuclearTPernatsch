Info << "Accelerating neutronic solution: implicit predictor (not tested for SP3)" << endl ;
{

    PtrList<volScalarField> leakages(energyGroups);
    forAll(leakages,energyI)
    {
        Dalbedo = D[energyI];
        fluxStarAlbedo = fluxStar2[energyI]*(SP3Neutronics);
        leakages.set
        (
            energyI,
            new volScalarField
            (
                IOobject
                (
                    "leakages"+Foam::name(energyI),
                    runTime.timeName(),
                    neutroMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::laplacian(D[energyI]/discFactor[energyI],fluxStar[energyI])/fluxStar[energyI]
            )
        );
    }

    PtrList<volScalarField> ddtFlux2(energyGroups);
    
    forAll(ddtFlux2,energyI)
    {
        ddtFlux2.set
        (
            energyI,
            new volScalarField
            (
                IOobject
                (
                    "ddtFlux2"+Foam::name(energyI),
                    runTime.timeName(),
                    neutroMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::ddt(IV[energyI]/discFactor[energyI],fluxStar2[energyI])
            )
        );
    }


    forAll(fluxStar[0].internalField(),i)
    {
        SquareMatrix<scalar> matrix(energyGroups,energyGroups,Foam::zero());
        List<scalar> source(energyGroups,0.0);
        List<scalar> predictedFlux(energyGroups,0.0);

        forAll(flux,energyI)//initialized predicted flux
        {
            predictedFlux[energyI] = fluxStar[energyI].internalField()[i] ;
        }   
        
        forAll(flux,energyI)
        {
            forAll(flux,energyJ)
            {
                if(energyJ==energyI)
                {
                    matrix[energyI][energyJ]=
                                  IV[energyI].internalField()[i]/runTime.deltaT().value()/discFactor[energyI].internalField()[i]
                                - nuSigmaEff[energyI].internalField()[i]/keff*(1.0-BetaTot.internalField()[i])*chiPrompt[energyI].internalField()[i]/discFactor[energyI].internalField()[i]
                                + sigmaDisapp[energyI].internalField()[i]/discFactor[energyI].internalField()[i]
                                - leakages[energyI].internalField()[i];
                }else
                {
                    matrix[energyI][energyJ]= 
                                - sigmaFromTo[energyJ][energyI].internalField()[i]/discFactor[energyJ].internalField()[i]
                                - nuSigmaEff[energyJ].internalField()[i]/keff*(1.0-BetaTot.internalField()[i])*chiPrompt[energyI].internalField()[i]/discFactor[energyJ].internalField()[i];
                }
            }

            source[energyI]= 
                      IV[energyI].internalField()[i]*fluxStar[energyI].oldTime().internalField()[i]/runTime.deltaT().value()/discFactor[energyI].internalField()[i]
                    //+ leakages[energyI].internalField()[i]
                    + delayedNeutroSource.internalField()[i]*chiDelayed[energyI].internalField()[i]
                    + 2*sigmaDisapp[energyI].internalField()[i]*fluxStar2[energyI].internalField()[i]/discFactor[energyI].internalField()[i]*SP3Neutronics 
                    + 2*ddtFlux2[energyI].internalField()[i]*SP3Neutronics;
        }

        solve(predictedFlux,matrix,source);

        forAll(fluxStar,energyI)
        {

            fluxStar[energyI].primitiveFieldRef()[i] = predictedFlux[energyI];
        
        }   

    }

    forAll(flux,energyI)
    {
        Dalbedo = D[energyI];
        fluxStarAlbedo = fluxStar2[energyI]*(SP3Neutronics);
        fluxStar[energyI].correctBoundaryConditions();
        flux[energyI]=(fluxStar[energyI]-fluxStar2[energyI]*2.0*SP3Neutronics)/discFactor[energyI];
        flux[energyI].correctBoundaryConditions();
    
    }

    #include "initializeNeutroSource.H"
}


