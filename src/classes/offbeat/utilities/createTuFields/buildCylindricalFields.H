Info << nl << "Building cylindrical volSymmTensor fields:" << nl << endl;

// Initialize ptrList of cylindrical fields
PtrList<volSymmTensorField> cylFields;

// Recompose cylindrical components fields into a cyl volTensorFields
forAll(rFields, i)
{
    // Get volSymmTensorField name
    word vstfCylName(rFields[i].name());
    vstfCylName.replace("_radial","Cyl");

    // Create volSymmTensorField
    volSymmTensorField vstfCyl
    (
        IOobject
        (
            vstfCylName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("", rFields[i].dimensions() ,symmTensor::zero),
        "zeroGradient"
    );

    // Set radial, theta and axial component
    vstfCyl  = rFields[i]        * symmTensor(1,0,0,0,0,0);
    vstfCyl += thetaFields[i]    * symmTensor(0,0,0,1,0,0);
    vstfCyl += zFields[i]        * symmTensor(0,0,0,0,0,1);

    // Add this vstf to ptrlist
    cylFields.append(&vstfCyl);

    vstfCyl.correctBoundaryConditions();

    // write the field
    vstfCyl.write();

    Info << tab << vstfCylName << " field recomposed from its components."
            << endl;
}