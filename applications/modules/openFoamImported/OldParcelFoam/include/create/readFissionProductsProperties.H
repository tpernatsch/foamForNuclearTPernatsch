// Populates the fission product species constants (FPYield_, FPLambda_,
// FPDp_, FPRho_, FPPd_) from constant/nuclearProperties. Included only when
// fpTransport_ is true, i.e. when a "fissionProducts" subDict is present.

{
    Info<< "\nReading fission products constants" << endl;

    const dictionary& fissionProductsDict = nuclearPropertiesDict_.subDict("fissionProducts");

    nFPSpecies_ = readInt(fissionProductsDict.lookup("species"));
    FPYield_.setSize(nFPSpecies_);
    FPLambda_.setSize(nFPSpecies_);
    FPDp_.setSize(nFPSpecies_);
    FPRho_.setSize(nFPSpecies_);
    FPPd_.setSize(nFPSpecies_);

    List<scalar> FPYieldRaw(fissionProductsDict.lookup("yield"));
    List<scalar> FPLambdaRaw(fissionProductsDict.lookup("lambda"));
    List<scalar> FPDpRaw(fissionProductsDict.lookup("dp"));
    List<scalar> FPRhoRaw(fissionProductsDict.lookup("rho"));
    List<scalar> FPPdRaw(fissionProductsDict.lookup("Pd"));

    if (FPLambdaRaw.size() == nFPSpecies_ && FPYieldRaw.size() == nFPSpecies_)
    {
        Info<< "    Found " << nFPSpecies_ << " species of fission products" << endl;
    }
    else
    {
        FatalError
            << "in " << fissionProductsDict.name()
            << " , invalid number of fission product species constants" << exit(FatalError);
    }

    forAll(FPYieldRaw, i)
    {
        FPYield_.set(i, new dimensionedScalar("FPYield" + Foam::name(i + 1), dimless, FPYieldRaw[i]));
        FPLambda_.set(i, new dimensionedScalar("FPLambda" + Foam::name(i + 1), dimless/dimTime, FPLambdaRaw[i]));
        FPDp_.set(i, new dimensionedScalar("FPDp" + Foam::name(i + 1), dimLength, FPDpRaw[i]));
        FPRho_.set(i, new dimensionedScalar("FPRho" + Foam::name(i + 1), dimDensity, FPRhoRaw[i]));
        FPPd_.set(i, new dimensionedScalar("FPPd" + Foam::name(i + 1), dimEnergy, FPPdRaw[i]));
    }
}
