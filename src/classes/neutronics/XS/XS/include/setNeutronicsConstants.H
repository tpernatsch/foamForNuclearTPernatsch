// This file sets constant values, that is, values that are not affected
// by feedbacks or at least are assumed to not change with time during
// the simulation.

Info<< "Set constant nuclear data (1/V, chiPrompt, chiDelayed, betaTot)"
    << endl;

forAll(nuSigmaEff_, energyI)
{
    Info<< "  Energy group " << energyI << endl;

    forAll(referenceZones_, zoneI)
    {
        label zone = zoneI;
        const word& name = referenceZones_[zoneI].keyword();
        label zoneId = mesh_.cellZones().findZoneID(name);

        Info<< "    Mesh zone: " << name << endl;

        if (zoneId == -1)
        {
            Warning
                << "mesh zone: " << name << " does not exists" << nl
                << "Ignore this warning if it is a completely withdrawn control rod"
                << endl;
        }
        else
        {
            forAll(mesh_.cellZones()[zoneId], cellIlocal)
            {
                label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
                zone = zoneI;

                // Constant
                IV_[energyI][cellIglobal] = IVList_[zone][energyI];
                // chiPrompt_[energyI][cellIglobal] = chiPromptList_[zone][energyI];
                // chiDelayed_[energyI][cellIglobal] = chiDelayedList_[zone][energyI];
                BetaTot_[cellIglobal] = BetaTotList_[zone];

                fuelFraction_[cellIglobal] = fuelFractionList_[zone];
                secondaryPowerVolumeFraction_[cellIglobal] = secondaryPowerVolumeFractionList_[zone];
                fractionToSecondaryPower_[cellIglobal] = fractionToSecondaryPowerList_[zone];

                // Variables
                D_[energyI][cellIglobal] = DList_[zone][energyI].getRef();
                nuSigmaEff_[energyI][cellIglobal] = nuSigmaEffList_[zone][energyI].getRef();
                sigmaPow_[energyI][cellIglobal] = sigmaPowList_[zone][energyI].getRef();
                sigmaRemoval_[energyI][cellIglobal] = sigmaRemovalList_[zone][energyI].getRef();
                chiPrompt_[energyI][cellIglobal] = chiPromptList_[zone][energyI].getRef();
                chiDelayed_[energyI][cellIglobal] = chiDelayedList_[zone][energyI].getRef();

                forAll(sigmaFromTo_, momentI)
                {
                    forAll(sigmaFromTo_[momentI], energyJ)
                    {
                        sigmaFromTo_[momentI][energyJ][energyI][cellIglobal] =
                            sigmaFromToList_[zone][momentI][energyJ][energyI].getRef();
                    }
                }
            }
        }
    }

    // Correct boundary conditions
    IV_[energyI].correctBoundaryConditions();
    BetaTot_.correctBoundaryConditions();
    fuelFraction_.correctBoundaryConditions();
    secondaryPowerVolumeFraction_.correctBoundaryConditions();
    fractionToSecondaryPower_.correctBoundaryConditions();
}

// Check that XSs have been specified for all regions
{
    forAll(mesh_.cellZones().names(), cellZoneI)
    {
        bool zoneInEntries(false);
        word zoneName(mesh_.cellZones().names()[cellZoneI]);
        forAll(referenceZones_, zoneI)
        {
            if (zoneName == referenceZones_[zoneI].keyword())
            {
                zoneInEntries = true;
                break;
            }
        }
        if (!zoneInEntries)
        {
            FatalErrorInFunction
                << "Zone: " << zoneName << " not found in neutronicProperties!"
                << exit(FatalError);
        }
    }
}
