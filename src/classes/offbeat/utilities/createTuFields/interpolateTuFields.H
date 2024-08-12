// Initialize list of Tu fields component of tensor field (i.e. sigma_r)
PtrList<volScalarField> rFields;
PtrList<volScalarField> thetaFields;
PtrList<volScalarField> zFields;

Info << "Translating fields Tu -> OFFBEAT:" << nl << endl;
forAll( TuFieldNames, i)
{
    if( not Tu2OFFBEATtable.found(TuFieldNames[i]) )
    {
        FatalErrorInFunction()
        << "Field name '" << TuFieldNames[i] << "' not found in Tu2OFFBEAT table!"
        << nl
        << " Please add it in Tu2OFFBEATtable.H HashTable if you want it to be "
        << " converted in an OFFBEAT field." << nl
        << abort(FatalError);
    }

    TuField& tufieldI = Tu2OFFBEATtable.find(TuFieldNames[i])(); 

    if ( tufieldI.nameOFFBEAT_ == "")
    {
        notTranslatedFields.append(TuFieldNames[i]);
    }
    else
    {
        // Get field names
        word tuFieldName(TuFieldNames[i]);
        word OFFBEATFieldName(tufieldI.nameOFFBEAT_);

        // Get unit conversion multiplier
        scalar uMultiplier(tufieldI.unitMultiplier_);

        // Print info
        Info << "  '"   << tuFieldName          << "'" << tab << " -> " 
            << "  '"    << OFFBEATFieldName     << "'"  
            << ", with conversion factor: "      << uMultiplier << "." 
            << endl;

        // Initialize volScalarField
        volScalarField vsf
        (
            IOobject
            (
                OFFBEATFieldName,
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("", tufieldI.dimension_, 0.0),
            "zeroGradient"
        );

        // Read Tu 2D map
        FieldField<Field, scalar> TuFieldData
        (
            PtrList<scalarField>
            (
                fieldsDict.lookup(tuFieldName),
                scalarFieldFieldINew()
            )
        );

        // Create 2D interpolation table
        scalarFieldInterpolateTable fieldInterpolationTable
        (
            zValues, 
            TuFieldData,
            Foam::interpolateTableBase::interpolationMethodNames_["linear"],
            interpolateTableBase::outofBoundsMehtodNames_["fixed"]
        );

        // Map the field on the the OFFBEAT mesh based on the interpolation tab
        forAll(mesh.C(), j)
        {
            // Get coordinates of this cell
            scalar cx(mesh.C()[j].x());
            scalar cy(mesh.C()[j].y());
            scalar cz(mesh.C()[j].z());

            // Compute the radius, assuming rod along z
            scalar radius(Foam::sqrt(cx*cx+cy*cy));

            // Create a 1D scalarField corrisponding to all the radial values
            // of the 2D field for this specific z
            scalarField rField(fieldInterpolationTable(cz));

            // Create a 1D interpolation table along r
            scalarInterpolateTable rTable
            (
                rValues, 
                rField, 
                interpolateTableBase::interpolationMethodNames_["linear"],
                interpolateTableBase::outofBoundsMehtodNames_["fixed"]
            );

            // Assign the interpolated value to this cell
            vsf[j] = rTable(radius);
        }

        // Use the multiplier to convert the field to OFFBEAT units
        vsf *= uMultiplier;

        // Special treatment for temperature conversion
        if (vsf.name() == "T")
        {
            vsf.field() += 273.15;
        }

        vsf.correctBoundaryConditions();

        // Check if the field is a radial, theta or axial component of a
        // volTensorField
        if
        ( 
            OFFBEATFieldName.find("_radial") != std::string::npos 
        )
        {
            rFields.append(&vsf);
        } 
        else if ( OFFBEATFieldName.find("_theta") != std::string::npos  )
        {
            thetaFields.append(&vsf);
        }
        else if ( OFFBEATFieldName.find("_axial") != std::string::npos )
        {
            zFields.append(&vsf);
        }
        else
        {
            // write the field in the apposite time folder
            vsf.write();
        }
    }
}