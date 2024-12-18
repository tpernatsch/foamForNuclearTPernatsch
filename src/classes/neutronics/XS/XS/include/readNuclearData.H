Info<< nl << "Reading nuclear data" << nl << endl;

// Log error if found previous GeN-Foam nuclearData file format
if (nuclearData_.found("zones"))
{
    FatalErrorInFunction()
        << "Found 'zones' dictionary in nuclearData file." << nl
        << "GeN-Foam now relies on the 'states' dictionary." << nl << nl
        << "states" << nl
        << "(" << nl
        << "    reference // Mandatory reference state" << nl
        << "    {" << nl
        << "        Tfuel    1000;" << nl << nl
        << "        zones" << nl
        << "        (" << nl
        << "            zone0" << nl
        << "            {" << nl
        << "                fuelFraction   ...;" << nl
        << "                ..." << nl
        << "            }" << nl
        << "            ..." << nl
        << "        );" << nl
        << "    }" << nl
        << "    perturbedState1" << nl
        << "    ..." << nl
        << ");" << nl << nl
        << "See examples in the GeN-Foam Tutorials" << nl
        << abort(FatalError);
}

// First XS state must be "reference"
if (states_.first().keyword() != "reference")
{
    FatalErrorInFunction()
        << "First XS state must be 'reference'. Found '"
        << states_.first().keyword() << "' as first state."
        << abort(FatalError);
}

Info<< nl << "Build XSs using "
    << energyGroups_ << " energy groups and "
    << precGroups_ << " precursor groups"
    << nl
    << endl;

Info<< "XS variables:" << endl;
scalarList variableRefValues(0);
forAll(xsVariablesDict_.toc(), varI)
{
    word variableName(xsVariablesDict_.toc()[varI]);

    word variableType(xsVariablesDict_.get<word>(variableName));

    xsVariableNames_[varI] = variableName;
    xsVariableTypes_[varI] = variableType;

    variableRefValues.append
    (
        referenceState_.get<scalar>(variableName)
    );

    Info<< tab
        << variableName
        << ", law: " << variableType
        << ", ref value: " << variableRefValues[varI]
        << endl;
}

// Reference perturbed parameters
// const scalar TfuelRef(referenceState_.lookupOrDefault("Tfuel", 0.01));
// const scalar TcladRef(referenceState_.lookupOrDefault("Tclad", 0.0));
// const scalar TcoolRef(referenceState_.lookupOrDefault("Tcool", 0.0));
// const scalar TstructMechRef(referenceState_.lookupOrDefault("TstructMech", 0.0));
// const scalar rhoCoolRef(referenceState_.lookupOrDefault("rhoCool", 0.0));
// const scalar axExpRef(referenceState_.lookupOrDefault("axExp", 0.0));
// const scalar radExpRef(referenceState_.lookupOrDefault("radExp", 0.0));

Info<< nl
    << "Use polyharmonic spline function mode " << polyharmonicSplineMode_
    << nl
    << endl;

// Look for reference state and extract reference values
Info<< "Zones:" << endl;
forAll(referenceZones_, zoneI)
{
    dictionary& dict = referenceZones_[zoneI].dict();
    const word& nameZone = referenceZones_[zoneI].keyword();

    Info<< tab << "Found zone: " << nameZone << endl;

    // Extract constant factors
    fuelFractionList_.set(zoneI, new scalar(dict.lookupOrDefault("fuelFraction", 1.0)));
    secondaryPowerVolumeFractionList_.set(zoneI, new scalar(dict.lookupOrDefault("secondaryPowerVolumeFraction", 1.0)));
    fractionToSecondaryPowerList_.set(zoneI, new scalar(dict.lookupOrDefault("fractionToSecondaryPower", 0.0)));
    dfAdjustList_.set(zoneI, new bool(dict.lookupOrDefault("dfAdjust", true)));

    // Nuclear data (assumed to be) not affected by feedbacks
    IVList_.set(zoneI, new scalarField("IV", dict, energyGroups_));
    chiPromptList_.set(zoneI, new scalarField("chiPrompt", dict, energyGroups_));
    chiDelayedList_.set(zoneI, new scalarField("chiDelayed", dict, energyGroups_));
    discFactorList_.set(zoneI, new scalarField("discFactor", dict, energyGroups_));
    integralFluxList_.set(zoneI, new scalarField("integralFlux", dict, energyGroups_));
    BetaList_.set(zoneI, new scalarField("Beta", dict, precGroups_));
    BetaTotList_.set(zoneI, new scalar(sum(BetaList_[zoneI])));
    lambdaList_.set(zoneI, new scalarField("lambda", dict, precGroups_));

    // Perturbed nuclear data
    DList_[zoneI].setSize(energyGroups_);
    nuSigmaEffList_[zoneI].setSize(energyGroups_);
    sigmaPowList_[zoneI].setSize(energyGroups_);
    sigmaRemovalList_[zoneI].setSize(energyGroups_);
    sigmaFromToList_[zoneI].setSize(legendreMoments_);

    forAll(DList_[zoneI], energyI)
    {
        DList_[zoneI].set
        (
            energyI,
            new nuclearDataOneEnergy(nameZone, polyharmonicSplineMode_, nxsVariables_)
        );
        nuSigmaEffList_[zoneI].set
        (
            energyI,
            new nuclearDataOneEnergy(nameZone, polyharmonicSplineMode_, nxsVariables_)
        );
        sigmaPowList_[zoneI].set
        (
            energyI,
            new nuclearDataOneEnergy(nameZone, polyharmonicSplineMode_, nxsVariables_)
        );
        sigmaRemovalList_[zoneI].set
        (
            energyI,
            new nuclearDataOneEnergy(nameZone, polyharmonicSplineMode_, nxsVariables_)
        );
    }

    forAll(sigmaFromToList_[zoneI], momentI)
    {
        sigmaFromToList_[zoneI][momentI].setSize(energyGroups_);

        forAll(sigmaFromToList_[zoneI][momentI], energyJ)
        {
            sigmaFromToList_[zoneI][momentI][energyJ].setSize(energyGroups_);

            forAll(sigmaFromToList_[zoneI][momentI][energyJ], energyI)
            {
                sigmaFromToList_[zoneI][momentI][energyJ].set
                (
                    energyI,
                    new nuclearDataOneEnergy(nameZone, polyharmonicSplineMode_, nxsVariables_)
                );
            }
        }
    }
}


// Read data
// Info<< "Parameters (Tfuel, Tclad, Tcool, TstructMech, rhoCool, axExp, radExp) << endl;"

forAll(states_, stateI)
{
    dictionary& dictState = states_[stateI].dict();
    const word& nameState = states_[stateI].keyword();

    // Extract perturbation input
    Info<< "State: " << nameState << " (";
    scalarList variableValues(0);
    forAll(xsVariableNames_, varI)
    {
        word variableName(xsVariableNames_[varI]);
        word variableType(xsVariableTypes_[varI]);

        scalar variableValue(dictState.lookupOrDefault(variableName, variableRefValues[varI]));

        // Print before transformation
        Info<< variableValue;
        if (varI < nxsVariables_-1)
        {
            Info<< ", ";
        }

        // Transformation
        if (variableType == "log")
        {
            variableValue = log(variableValue);
        }
        else if (variableType == "sqrt")
        {
            variableValue = sqrt(variableValue);
        }

        variableValues.append(variableValue);
    }
    Info<< ")" << endl;

    // const scalar Tfuel(dictState.lookupOrDefault("Tfuel", TfuelRef));
    // const scalar Tclad(dictState.lookupOrDefault("Tclad", TcladRef));
    // const scalar Tcool(dictState.lookupOrDefault("Tcool", TcoolRef));
    // const scalar TstructMech(dictState.lookupOrDefault("TstructMech", TstructMechRef));
    // const scalar rhoCool(dictState.lookupOrDefault("rhoCool", rhoCoolRef));
    // const scalar axExp(dictState.lookupOrDefault("axExp", axExpRef));
    // const scalar radExp(dictState.lookupOrDefault("radExp", radExpRef));

    // const scalar TfuelSpectrum(fastNeutrons_ ? log(Tfuel) : pow(Tfuel, 0.5));

    // Info<< "State: " << nameState
    //     << " ("
    //     << Tfuel << ", " << Tclad << ", " << Tcool << ", "
    //     << TstructMech << ", " << rhoCool << ", " << axExp << ", "
    //     << radExp
    //     << ")"
    //     << endl;

    // Extract zones from perturbed state
    PtrList<entry> zones = dictState.lookup("zones");

    forAll(zones, zoneI)
    {
        dictionary& dict = zones[zoneI].dict();
        const word& nameZone = zones[zoneI].keyword();

        // Look for the zone index in the reference zone list
        // The perturbed state zones list may have different length compared to
        // the reference zones list
        label zoneRefI(-1);
        forAll(referenceZones_, zoneII)
        {
            if (referenceZones_[zoneII].keyword() == nameZone)
            {
                zoneRefI = zoneII;
                break;
            }
        }

        if (zoneRefI == -1)
        {
            Warning << nl
                << "Zone name " << nameZone << " not found in reference state." << nl
                << "No perturbation will be applied for this zone." << nl
                << "Include this zone name in the 'reference' state or remove it" << nl
                << "from " << nameState << " state to remove this message" << nl
                << endl;
            continue;
        }

        Info<< tab << "zone name: " << nameZone << endl;

        const scalarField DList("D", dict, energyGroups_);
        const scalarField nuSigmaEffList("nuSigmaEff", dict, energyGroups_);
        const scalarField sigmaPowList("sigmaPow", dict, energyGroups_);
        const scalarField sigmaRemovalList("sigmaRemoval", dict, energyGroups_);


        // Fill interpolation objects per zones and per energy for all states
        forAll(DList_[zoneRefI], energyI)
        {
            DList_[zoneRefI][energyI].addData
            (
                DList[energyI],
                variableValues
                // {
                //     TfuelSpectrum, Tclad, Tcool, TstructMech,
                //     rhoCool, axExp, radExp
                // }
            );

            nuSigmaEffList_[zoneRefI][energyI].addData
            (
                nuSigmaEffList[energyI],
                variableValues
                // {
                //     TfuelSpectrum, Tclad, Tcool, TstructMech,
                //     rhoCool, axExp, radExp
                // }
            );

            sigmaPowList_[zoneRefI][energyI].addData
            (
                sigmaPowList[energyI],
                variableValues
                // {
                //     TfuelSpectrum, Tclad, Tcool, TstructMech,
                //     rhoCool, axExp, radExp
                // }
            );

            sigmaRemovalList_[zoneRefI][energyI].addData
            (
                sigmaRemovalList[energyI],
                variableValues
                // {
                //     TfuelSpectrum, Tclad, Tcool, TstructMech,
                //     rhoCool, axExp, radExp
                // }
            );
        }

        forAll(sigmaFromToList_[zoneRefI], momentI)
        {
            const scalarSquareMatrix sigmaFromToList
            (
                dict.lookup("scatteringMatrixP" + Foam::name(momentI))
            );

            forAll(DList_[zoneRefI], energyJ)
            {
                forAll(DList_[zoneRefI], energyI)
                {
                    sigmaFromToList_[zoneRefI][momentI][energyJ][energyI].addData
                    (
                        sigmaFromToList(energyJ, energyI),
                        variableValues
                        // {
                        //     TfuelSpectrum, Tclad, Tcool, TstructMech,
                        //     rhoCool, axExp, radExp
                        // }
                    );
                }
            }
        }
    }
}


Info<< "Number of control rods: " << CRNumber_ << endl;
forAll(CRentries_, CRzoneI)
{
    dictionary& CRdict = CRentries_[CRzoneI].dict();

    CRstart_.set(CRzoneI, new scalar(CRdict.lookupOrDefault("startTime", 0.0)));
    CRfinish_.set(CRzoneI, new scalar(CRdict.lookupOrDefault("endTime", 0.0)));
    CRspeed_.set(CRzoneI, new scalar(CRdict.lookupOrDefault("speed", 0.0)));
    CRFollowerName_.set(CRzoneI, new word(CRdict.lookup("followerName")));
    initialDistanceFromMeshCR_.set(CRzoneI, new scalar(CRdict.lookupOrDefault("initialDistanceFromMeshCR", 0.0)));

    label FOLLOWERzoneId = mesh.cellZones().findZoneID(CRFollowerName_[CRzoneI]);

    if (FOLLOWERzoneId == -1)
    {
        Warning
            << "Control rod with follower: " << CRFollowerName_[CRzoneI]
            << " does not exists"
            << endl;
    }

    scalar cellHeight(-100.0);
    scalar maxFollowerHeight(cellHeight);

    volScalarField centersVSF(mesh.C() & vector(0, 0, 1));

    Field<scalar> centersSF(centersVSF.internalField());

    forAll(mesh.cellZones()[FOLLOWERzoneId], cellIlocal)
    {
        label cellIglobal = mesh.cellZones()[FOLLOWERzoneId][cellIlocal];

        maxFollowerHeight = max(cellHeight, centersSF[cellIglobal]);

        cellHeight = maxFollowerHeight;
    }

    reduce(maxFollowerHeight, maxOp<scalar>());
    CRinitialPosition_.set(CRzoneI, new scalar(maxFollowerHeight - initialDistanceFromMeshCR_[CRzoneI]));

    if (mesh.time().value() < CRfinish_[CRzoneI])
    {
        CRinitialPosition_[CRzoneI] = CRinitialPosition_[CRzoneI]
            - CRspeed_[CRzoneI] * (max(mesh.time().value(), CRstart_[CRzoneI]) - CRstart_[CRzoneI]);
    }
    else
    {
        CRinitialPosition_[CRzoneI] = CRinitialPosition_[CRzoneI]
            - CRspeed_[CRzoneI] * (CRfinish_[CRzoneI] - CRstart_[CRzoneI]);
    }
    CRposition_.set(CRzoneI, new scalar(CRinitialPosition_[CRzoneI]));

    Info<< "Control rod " << CRentries_[CRzoneI].keyword()
        << " in position: " << CRposition_[CRzoneI]
        << endl;
}

if (isLowMemory_)
{
    forAll(nuSigmaEff_, energyI)
    {
        forAll(nuSigmaEff_, energyJ)
        {
            label exchangeYesNo(0);
            forAll(referenceZones_, zoneI)
            {
                if (sigmaFromToList_[zoneI][0][energyJ][energyI].getRef() > 5.0e-6)
                {
                    exchangeYesNo = 1;
                    break;
                }
            }
            sigmaFromToYesNo_[energyJ][energyI] = exchangeYesNo;
        }
    }
}


Info<< "Finished reading nuclear data " << endl;
