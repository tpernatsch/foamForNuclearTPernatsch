// Calculate Control Rod position
forAll(CRentries_, CRzoneI)
{
    // Recalculate positions in case mesh has moved
    label FOLLOWERzoneId = mesh_.cellZones().findZoneID(CRFollowerName_[CRzoneI]);

    if (FOLLOWERzoneId == -1)
    {
        Warning
            << "Control rod with follower: " << CRFollowerName_[CRzoneI]
            << " does not exists"
            << endl;
    }

    scalar cellHeight(-100.0);
    scalar maxFollowerHeight(cellHeight);

    volScalarField centersVSF(mesh_.C() & vector(0, 0, 1));

    Field<scalar> centersSF(centersVSF.internalField());

    forAll(mesh_.cellZones()[FOLLOWERzoneId], cellIlocal)
    {
        label cellIglobal = mesh_.cellZones()[FOLLOWERzoneId][cellIlocal];

        maxFollowerHeight = max(cellHeight, centersSF[cellIglobal]);

        cellHeight = maxFollowerHeight;
    }

    reduce(maxFollowerHeight, maxOp<scalar>());
    CRinitialPosition_[CRzoneI] = maxFollowerHeight - initialDistanceFromMeshCR_[CRzoneI];

    if (mesh_.time().value() < CRfinish_[CRzoneI])
    {
        CRposition_[CRzoneI] = CRinitialPosition_[CRzoneI]
            - CRspeed_[CRzoneI] * (max(mesh_.time().value(), CRstart_[CRzoneI]) - CRstart_[CRzoneI]);
    }
    else
    {
        CRposition_[CRzoneI] = CRinitialPosition_[CRzoneI]
            - CRspeed_[CRzoneI] * (CRfinish_[CRzoneI] - CRstart_[CRzoneI]);
    }

    Info<< "Control rod " << CRentries_[CRzoneI].keyword()
        << " in position: " << CRposition_[CRzoneI]
        << endl;
}


// // Compute expansions
// tmp<volVectorField> tDNeutroRad(new volVectorField(Disp));
// volVectorField& DNeutroRad = tDNeutroRad.ref();

// tmp<volScalarField> tDNeutroAx(new volScalarField(Disp & axialOrientation_));
// volScalarField& DNeutroAx = tDNeutroAx.ref();

// DNeutroRad.replace(0, DNeutroRad.component(0)*(1.0 - axialOrientation_[0]));
// DNeutroRad.replace(1, DNeutroRad.component(1)*(1.0 - axialOrientation_[1]));
// DNeutroRad.replace(2, DNeutroRad.component(2)*(1.0 - axialOrientation_[2]));

// DNeutroRad.correctBoundaryConditions();
// DNeutroAx.correctBoundaryConditions();

// radExp_ = fvc::div(DNeutroRad)/2.0;
// axExp_ = fvc::grad(DNeutroAx) & axialOrientation_;

// radExp_.correctBoundaryConditions();
// axExp_.correctBoundaryConditions();

// tDNeutroRad.clear();
// tDNeutroAx.clear();

// log(T) if fast spectrum or sqrt(T) if thermal
// logT_ = log(max(0.01, Tfuel/dimensionedScalar("", dimTemperature, 1.0)));
// sqrtT_ = pow(max(dimensionedScalar("", dimTemperature, 0.0), Tfuel), 0.5);

// Change variable law
forAll(xsVariableNames_, varI)
{
    word variableName(xsVariableNames_[varI]);
    word variableType(xsVariableTypes_[varI]);

    // Compute expansion
    if (variableName == "radExp")
    {
        tmp<volVectorField> tDNeutroRad(new volVectorField(disp_));
        volVectorField& DNeutroRad = tDNeutroRad.ref();

        DNeutroRad.replace(0, DNeutroRad.component(0)*(1.0 - axialOrientation_[0]));
        DNeutroRad.replace(1, DNeutroRad.component(1)*(1.0 - axialOrientation_[1]));
        DNeutroRad.replace(2, DNeutroRad.component(2)*(1.0 - axialOrientation_[2]));

        DNeutroRad.correctBoundaryConditions();

        // Compute radial expansion
        xsVariableFields_[varI] = fvc::div(DNeutroRad)/2.0;

        xsVariableFields_[varI].correctBoundaryConditions();

        tDNeutroRad.clear();
    }
    else if (variableName == "axExp")
    {
        tmp<volScalarField> tDNeutroAx(new volScalarField(disp_ & axialOrientation_));
        volScalarField& DNeutroAx = tDNeutroAx.ref();

        DNeutroAx.correctBoundaryConditions();

        // Compute axial expansion
        xsVariableFields_[varI] = fvc::grad(DNeutroAx) & axialOrientation_;

        xsVariableFields_[varI].correctBoundaryConditions();

        tDNeutroAx.clear();
    }

    // Tranform with user-specified law
    if (variableType == "log")
    {
        xsVariableTransformedFields_.set(varI, log
        (
            max
            (
                xsVariableFields_[varI] / dimensionedScalar("", xsVariableFields_[varI].dimensions(), 1.0),
                0.01
            )
        ));
    }
    else if (variableType == "sqrt")
    {
        xsVariableTransformedFields_.set(varI, pow
        (
            max
            (
                xsVariableFields_[varI],
                dimensionedScalar("", xsVariableFields_[varI].dimensions(), 0.0)
            ),
            0.5
        ));
    }
}


scalarList variableValues(nxsVariables_);

forAll(nuSigmaEff_, energyI)
{
    bool isParametrize(true);
    forAll(doNotParametrize_, groupI)
    {
        if (energyI == doNotParametrize_[groupI])
        {
            isParametrize = false;
            Info<< "Not Parametrizing XS for group: " << doNotParametrize_[groupI]
                << endl;
            break;
        }
    }

    forAll(referenceZones_, zoneI)
    {
        label zone = zoneI;
        const word& name = referenceZones_[zoneI].keyword();
        label zoneId = mesh_.cellZones().findZoneID(name);

        if (zoneId == -1)
        {
            Warning
                << "Mesh zone " << name << " does not exists" << nl
                << "Ignore this warning if it is a completely withdrawn control rod"
                << endl;

            continue;
        }

        // Check if you are in a follower or CR
        bool isFollower(false);
        bool isCR(false);
        label CRzone(0);

        forAll(CRentries_, CRzoneI)
        {
            const word& CRname = CRentries_[CRzoneI].keyword();

            if (name == CRname)
            {
                isCR = true;
                CRzone = CRzoneI;
            }
            if (name == CRFollowerName_[CRzoneI])
            {
                isFollower = true;
                CRzone = CRzoneI;
            }
        }

        forAll(mesh_.cellZones()[zoneId], cellIlocal)
        {
            label cellIglobal = mesh_.cellZones()[zoneId][cellIlocal];
            zone = zoneI;

            if (CRNumber_ > 0)
            {
                if (isFollower && mesh_.C()[cellIglobal].z() > CRposition_[CRzone])
                {
                    forAll(referenceZones_, zoneII)
                    {
                        if (referenceZones_[zoneII].keyword() == CRentries_[CRzone].keyword())
                        {
                            zone = zoneII;
                            break;
                        }
                    }
                }

                if (isCR && mesh_.C()[cellIglobal].z() < CRposition_[CRzone])
                {
                    forAll(referenceZones_, zoneII)
                    {
                        if (referenceZones_[zoneII].keyword() == CRFollowerName_[CRzone])
                        {
                            zone = zoneII;
                            break;
                        }
                    }
                }
            }

            // Adjustment of D and XS for current reactor conditions of:
            // Tfuel, cladExp, Tcool, rhoCool, TStructMech, axialExp and
            // radialExp.

            // scalarList variableValues({
            //     fastNeutrons_? logT_[cellIglobal] : sqrtT_[cellIglobal],
            //     Tclad[cellIglobal],
            //     Tcool[cellIglobal],
            //     TStructMech[cellIglobal],
            //     rhoCool[cellIglobal],
            //     axExp_[cellIglobal],
            //     radExp_[cellIglobal]
            // });

            if (!isReadWriteXS_)
            {
                forAll(variableValues, varI)
                {
                    if (xsVariableTypes_[varI] == "lin")
                    {
                        variableValues[varI] = xsVariableFields_[varI][cellIglobal];
                    }
                    else
                    {
                        variableValues[varI] = xsVariableTransformedFields_[varI][cellIglobal];
                    }
                }

                D_[energyI][cellIglobal] = DList_[zone][energyI].get
                (
                    variableValues,
                    isParametrize
                );

                nuSigmaEff_[energyI][cellIglobal] = nuSigmaEffList_[zone][energyI].get
                (
                    variableValues,
                    isParametrize
                );

                sigmaPow_[energyI][cellIglobal] = sigmaPowList_[zone][energyI].get
                (
                    variableValues,
                    isParametrize
                );

                sigmaRemoval_[energyI][cellIglobal] = sigmaRemovalList_[zone][energyI].get
                (
                    variableValues,
                    isParametrize
                );

                if (!isLowMemory_)
                {
                    forAll(sigmaFromTo_, momentI)
                    {
                        forAll(sigmaFromTo_[momentI], energyJ)
                        {
                            sigmaFromTo_[momentI][energyJ][energyI][cellIglobal] =
                                sigmaFromToList_[zone][momentI][energyJ][energyI].get
                                (
                                    variableValues,
                                    isParametrize
                                );
                        }
                    }
                }
            }

            if (
                useGivenDiscFactors_ ||
                (mesh_.time().timeIndex() - mesh_.time().startTimeIndex() > 0)
            )
            {
                discFactor_[energyI][cellIglobal] = discFactorList_[zone][energyI];
            }
            else
            {
                discFactorList_[zone][energyI] = discFactor_[energyI][cellIglobal];
            }

        } // close forall cellIlocal in zoneId
    } // close forall zoneI in zone referenceZones_


    // Correct boundary conditions
    D_[energyI].correctBoundaryConditions();
    nuSigmaEff_[energyI].correctBoundaryConditions();
    sigmaPow_[energyI].correctBoundaryConditions();
    sigmaRemoval_[energyI].correctBoundaryConditions();
    forAll(sigmaFromTo_, momentI)
    {
        forAll(sigmaFromTo_[momentI], energyJ)
        {
            sigmaFromTo_[momentI][energyJ][energyI].correctBoundaryConditions();
        }
    }
    discFactor_[energyI].correctBoundaryConditions();
}
