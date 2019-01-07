
{

    Info << "Accelerating neutronic solution: integral predictor" << endl ;

    oneGroupFlux *= 0.0;
    forAll(fluxStar, energyI)
    {
           oneGroupFlux +=fluxStar[energyI];
    }
    dimensionedScalar integralFlux = fvc::domainIntegrate(oneGroupFlux);

    volScalarField oneGroupFluxOld("",oneGroupFlux);
    oneGroupFluxOld *= 0.0;
    forAll(fluxStar, energyI)
    {
           oneGroupFluxOld +=fluxStar[energyI].oldTime();
    }

    dimensionedScalar integralFluxOld = fvc::domainIntegrate(oneGroupFluxOld);

    dimensionedScalar temporalTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar temporalTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar leakageTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar leakageTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar disappearenceTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar disappearenceTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar sourceTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar sourceTerm("",(dimless/dimLength),scalar(0.0));

    dimensionedScalar scatteringSourceTermInt("",(dimless/dimTime),scalar(0.0));
    dimensionedScalar scatteringSourceTerm("",(dimless/dimLength),scalar(0.0));

    for(label energyI=0; energyI<energyGroups ; energyI++)
    {

        Dalbedo = D[energyI];

        temporalTermInt += fvc::domainIntegrate(IV[energyI]*fluxStar[energyI]/discFactor[energyI]/mesh_.time().deltaT());

        leakageTermInt += fvc::domainIntegrate(fvc::laplacian(D[energyI]/discFactor[energyI], fluxStar[energyI]));

        disappearenceTermInt += fvc::domainIntegrate(sigmaDisapp[energyI]*fluxStar[energyI]/discFactor[energyI]);

        #include "calcScatteringSource.H"
        scatteringSourceTermInt += fvc::domainIntegrate(scatteringSourceExtra);
    }

    sourceTermInt = fvc::domainIntegrate(neutroSource/keff_*(1.0-BetaTot));

    dimensionedScalar precursorTermInt = fvc::domainIntegrate(delayedNeutroSource);

    temporalTerm = temporalTermInt/integralFlux;

    leakageTerm = leakageTermInt/integralFlux;

    disappearenceTerm = disappearenceTermInt/integralFlux;

    sourceTerm = sourceTermInt/integralFlux;

    scatteringSourceTerm = scatteringSourceTermInt/integralFlux;

    dimensionedScalar precursorTerm = precursorTermInt;

    Info<< "Predict new integral flux " << endl;
    dimensionedScalar newIntegralFlux = (precursorTerm + integralFluxOld*temporalTerm)/(temporalTerm-leakageTerm+disappearenceTerm-sourceTerm-scatteringSourceTerm);

    Info<< "Old integral flux: " << integralFlux.value() << endl;
    Info<< "New predicted integral flux: " << newIntegralFlux.value() << endl;

    Info<< "Reconstruct fluxes " << endl;
    forAll(fluxStar, energyI)
    {
        Dalbedo = D[energyI];
        fluxStarAlbedo *= 0.0;
        fluxStar[energyI] *= newIntegralFlux/integralFluxOld;
        fluxStar[energyI].correctBoundaryConditions();
        flux[energyI]=fluxStar[energyI]/discFactor[energyI];
        flux[energyI].correctBoundaryConditions();
    }

    #include "initializeNeutroSource.H"
}

