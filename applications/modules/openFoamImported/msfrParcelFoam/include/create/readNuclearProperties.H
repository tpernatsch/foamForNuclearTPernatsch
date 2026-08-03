// Populates the group-indexed neutron transport, delayed-neutron precursor
// and decay-heat precursor constants (D_ref_, Sa_ref_, ..., precLambda_,
// decLambda_, ...) read from constant/nuclearProperties. Included once from
// the constructor body, after nEnergyGroups_/nPrecGroups_/nDecGroups_ have
// already been read (see the initializer list).

{
    Info<< "\nReading neutron transport constants" << endl;

    const dictionary& neutronTransportDict = nuclearPropertiesDict_.subDict("neutronTransport");

    List<scalar> D_refRaw(neutronTransportDict.lookup("DRef"));
    List<scalar> alphaD_refRaw(neutronTransportDict.lookup("alphaDRef"));

    List<scalar> Sa_refRaw(neutronTransportDict.lookup("SaRef"));
    List<scalar> alphaSa_refRaw(neutronTransportDict.lookup("alphaSaRef"));

    List<scalar> Sf_refRaw(neutronTransportDict.lookup("SfRef"));
    List<scalar> alphaSf_refRaw(neutronTransportDict.lookup("alphaSfRef"));

    List<scalar> Ss_refRaw(neutronTransportDict.lookup("SsRef"));
    List<scalar> alphaSs_refRaw(neutronTransportDict.lookup("alphaSsRef"));

    List<scalar> NuRaw(neutronTransportDict.lookup("Nu"));
    List<scalar> EfRaw(neutronTransportDict.lookup("Ef"));
    List<scalar> invVelRaw(neutronTransportDict.lookup("inverseVelocity"));
    List<scalar> chiPromptRaw(neutronTransportDict.lookup("promptYeld"));
    List<scalar> chiDelayedRaw(neutronTransportDict.lookup("delayedYeld"));

    if (D_refRaw.size() == nEnergyGroups_)
    {
        Info<< "    Found " << nEnergyGroups_ << " energy groups" << endl;
    }
    else
    {
        FatalError
            << "in " << neutronTransportDict.name()
            << " , invalid number of neutron transport constants" << exit(FatalError);
    }

    forAll(D_refRaw, i)
    {
        word iName = Foam::name(i + 1);

        D_ref_.set(i, new dimensionedScalar("D" + iName + "_ref", dimLength, D_refRaw[i]));
        Sa_ref_.set(i, new dimensionedScalar("Sa" + iName + "_ref", dimless/dimLength, Sa_refRaw[i]));
        Sf_ref_.set(i, new dimensionedScalar("Sf" + iName + "_ref", dimless/dimLength, Sf_refRaw[i]));

        alphaD_ref_.set(i, new dimensionedScalar("alphaD" + iName + "_ref", dimLength, alphaD_refRaw[i]));
        alphaSa_ref_.set(i, new dimensionedScalar("alphaSa" + iName + "_ref", dimless/dimLength, alphaSa_refRaw[i]));
        alphaSf_ref_.set(i, new dimensionedScalar("alphaSf" + iName + "_ref", dimless/dimLength, alphaSf_refRaw[i]));

        forAll(D_refRaw, j)
        {
            label ijInd = nEnergyGroups_*i + j;
            word ijName = Foam::name(i + 1) + Foam::name(j + 1);

            Ss_ref_.set(ijInd, new dimensionedScalar("Ss" + ijName + "_ref", dimless/dimLength, Ss_refRaw[ijInd]));
            alphaSs_ref_.set(ijInd, new dimensionedScalar("alphaSs" + ijName + "_ref", dimless/dimLength, alphaSs_refRaw[ijInd]));
        }

        Nu_.set(i, new dimensionedScalar("Nu" + Foam::name(i + 1), dimless, NuRaw[i]));
        Ef_.set(i, new dimensionedScalar("Ef" + Foam::name(i + 1), dimEnergy, EfRaw[i]));
        invVel_.set(i, new dimensionedScalar("invVel" + Foam::name(i + 1), dimTime/dimLength, invVelRaw[i]));
        chiPrompt_.set(i, new dimensionedScalar("chiPrompt" + Foam::name(i + 1), dimless, chiPromptRaw[i]));
        chiDelayed_.set(i, new dimensionedScalar("chiDelayed" + Foam::name(i + 1), dimless, chiDelayedRaw[i]));
    }

    Info<< "\nReading delayed neutron precursors constants" << endl;

    const dictionary& neutronPrecursorsDict = nuclearPropertiesDict_.subDict("delayedNeutronPrecursors");

    List<scalar> precLambdaRaw(neutronPrecursorsDict.lookup("lambda"));
    List<scalar> precBetaRaw(neutronPrecursorsDict.lookup("beta"));

    if (precLambdaRaw.size() == nPrecGroups_ && precBetaRaw.size() == nPrecGroups_)
    {
        Info<< "    Found " << nPrecGroups_ << " groups of delayed neutron precursors" << endl;
    }
    else
    {
        FatalError
            << "in " << neutronPrecursorsDict.name()
            << " , invalid number of delayed neutron precursors constants" << exit(FatalError);
    }

    forAll(precLambdaRaw, i)
    {
        precLambda_.set(i, new dimensionedScalar("precLambda" + Foam::name(i + 1), dimless/dimTime, precLambdaRaw[i]));
        precBeta_.set(i, new dimensionedScalar("precBeta" + Foam::name(i + 1), dimless, precBetaRaw[i]));
        precBetaTot_ += precBeta_[i];
    }

    Info<< "\nReading decay heat precursors constants" << endl;

    const dictionary& decayHeatPrecursorsDict = nuclearPropertiesDict_.subDict("decayHeatPrecursors");

    List<scalar> decLambdaRaw(decayHeatPrecursorsDict.lookup("lambda"));
    List<scalar> decBetaRaw(decayHeatPrecursorsDict.lookup("beta"));

    if (decLambdaRaw.size() == nDecGroups_ && decBetaRaw.size() == nDecGroups_)
    {
        Info<< "    Found " << nDecGroups_ << " groups of decay heat precursors" << endl;
    }
    else
    {
        FatalError
            << "in " << decayHeatPrecursorsDict.name()
            << " , invalid number of decay heat precursors constants" << exit(FatalError);
    }

    forAll(decLambdaRaw, i)
    {
        decLambda_.set(i, new dimensionedScalar("decLambda" + Foam::name(i + 1), dimless/dimTime, decLambdaRaw[i]));
        decBeta_.set(i, new dimensionedScalar("decBeta" + Foam::name(i + 1), dimless, decBetaRaw[i]));
        decBetaTot_ += decBeta_[i];
    }
}
