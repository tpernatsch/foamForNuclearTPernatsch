double legendrePoly(int, double);

Info<< "Calculating Legendre matrices" << nl << endl;

Info<< "Warning: considering Legendre scattering moments up to 7th at maximum." << endl;

scalarSquareMatrix muMatrix(discreteDirections_);

scalarSquareMatrix weightMatrix(discreteDirections_);

forAll(directionVersors_, directionI)
{
    forAll(directionVersors_, directionJ)
    {
        muMatrix[directionI][directionJ] = directionVersors_[directionI] & directionVersors_[directionJ];
        weightMatrix[directionI][directionJ] = directionWeights_[directionJ];
    }
}

Info<< "Done\n" << endl;

PtrList<scalar> maxLegendreIntegral(xs_.legendreMoments());
PtrList<scalar> minLegendreIntegral(xs_.legendreMoments());
PtrList<scalar> maxAbsLegendreIntegral(xs_.legendreMoments());
PtrList<scalar> minAbsLegendreIntegral(xs_.legendreMoments());
scalar legendreIntegral;
scalar legendreAbsIntegral;

forAll(legendreMatrices_, legendreI)
{
    Info<< "Calculating Legendre matrices... Legendre moment: " << legendreI << endl;
    maxLegendreIntegral.set(legendreI, new scalar(-100.0));
    minLegendreIntegral.set(legendreI, new scalar(100.0));
    maxAbsLegendreIntegral.set(legendreI, new scalar(-100.0));
    minAbsLegendreIntegral.set(legendreI, new scalar(100.0));

    // Info<< "Also, checking consistency of the matrix" << endl;
    legendreMatrices_.set
    (
        legendreI,
        new scalarSquareMatrix(discreteDirections_)
    );
    forAll(directionVersors_, directionI)
    {
        legendreIntegral = 0;
        legendreAbsIntegral = 0;
        forAll(directionVersors_, directionJ)
        {
            legendreMatrices_[legendreI][directionI][directionJ] = legendrePoly(legendreI, muMatrix[directionI][directionJ]);
            legendreIntegral += legendreMatrices_[legendreI][directionI][directionJ] * directionWeights_[directionJ] * Foam::constant::mathematical::pi * 4.0;
            legendreAbsIntegral += fabs(legendreMatrices_[legendreI][directionI][directionJ]) * directionWeights_[directionJ] * Foam::constant::mathematical::pi * 4.0;
        }
        maxLegendreIntegral[legendreI] = max(maxLegendreIntegral[legendreI],legendreIntegral);
        minLegendreIntegral[legendreI] = min(minLegendreIntegral[legendreI],legendreIntegral);
        maxAbsLegendreIntegral[legendreI] = max(maxAbsLegendreIntegral[legendreI],legendreAbsIntegral);
        minAbsLegendreIntegral[legendreI] = min(minAbsLegendreIntegral[legendreI],legendreAbsIntegral);
    }
}

Info<< "Done" << nl << endl;

Info<< "Consistency report on Legendre Integrals:" << nl
    << "Moment |   Minimum    |   Maximum    | Minimum Abs  | Maximum Abs  |    Theory   " << nl
    << "-------|--------------|--------------|--------------|--------------|-------------"
    << endl;

List<scalar> theoreticalLegendreAbsIntegral({12.5664, 6.28319, 4.8366, 4.08407, 3.6014, 3.25777, 2.997, 2.79031});

Info.stream().setf(ios_base::left);

forAll(legendreMatrices_, legendreI)
{
    Info<< "P" << setw(5) << legendreI << " | "
        << setw(12) << minLegendreIntegral[legendreI] << " | "
        << setw(12) << maxLegendreIntegral[legendreI] << " | "
        << setw(12) << minAbsLegendreIntegral[legendreI] << " | "
        << setw(12) << maxAbsLegendreIntegral[legendreI] << " | "
        << setw(12) << theoreticalLegendreAbsIntegral[legendreI]
        << endl;
}
