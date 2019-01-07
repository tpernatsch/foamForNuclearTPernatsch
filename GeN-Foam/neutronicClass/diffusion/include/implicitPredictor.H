Info << "Accelerating neutronic solution: implicit predictor (not tested for SP3)" << endl ;
{

    PtrList<volScalarField> leakages(xs.energyGroups());
    forAll(leakages,energyI)
    {
        Dalbedo = D[energyI];
        fluxStarAlbedo *= 0.0;
        leakages.set
        (
            energyI,
            new volScalarField
            (
                IOobject
                (
                    "leakages"+Foam::name(energyI),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::laplacian(D[energyI]/discFactor[energyI],fluxStar[energyI])/fluxStar[energyI]
            )
        );
    }

    PtrList<volScalarField> ddtFlux2(xs.energyGroups());


    forAll(fluxStar[0].internalField(),i)
    {
        SquareMatrix<scalar> matrix(xs.energyGroups(),xs.energyGroups(),Foam::zero());
        List<scalar> source(xs.energyGroups(),0.0);
        List<scalar> predictedFlux(xs.energyGroups(),0.0);

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
                                  IV[energyI].internalField()[i]/mesh_.time().deltaT().value()/discFactor[energyI].internalField()[i]
                                - nuSigmaEff[energyI].internalField()[i]/keff_*(1.0-BetaTot.internalField()[i])*chiPrompt[energyI].internalField()[i]/discFactor[energyI].internalField()[i]
                                + sigmaDisapp[energyI].internalField()[i]/discFactor[energyI].internalField()[i]
                                - leakages[energyI].internalField()[i];
                }else
                {
                    matrix[energyI][energyJ]= 
                                - sigmaFromTo[energyJ][energyI].internalField()[i]/discFactor[energyJ].internalField()[i]
                                - nuSigmaEff[energyJ].internalField()[i]/keff_*(1.0-BetaTot.internalField()[i])*chiPrompt[energyI].internalField()[i]/discFactor[energyJ].internalField()[i];
                }
            }

            source[energyI]= 
                      IV[energyI].internalField()[i]*fluxStar[energyI].oldTime().internalField()[i]/mesh_.time().deltaT().value()/discFactor[energyI].internalField()[i]
                    //+ leakages[energyI].internalField()[i]
                    + delayedNeutroSource.internalField()[i]*chiDelayed[energyI].internalField()[i];
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
        fluxStarAlbedo *= 0.0;
        fluxStar[energyI].correctBoundaryConditions();
        flux[energyI]=fluxStar[energyI]/discFactor[energyI];
        flux[energyI].correctBoundaryConditions();
    
    }

    #include "initializeNeutroSource.H"
}


