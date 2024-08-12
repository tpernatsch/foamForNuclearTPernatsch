Info << nl << "Building cartesian volSymmTensor fields:" << nl << endl;

// The carthesian fields are obtained through the transformation:
//
//             T_cartesian = J * T_cylindrical * J^T
//
// Where the jacobian is:
//      
//                  | cos(theta) -rsin(theta)  0 |
//              J = | sin(theta)  rcos(theta)  0 |
//                  | 0               0        1 |
//  
// Since cos(theta)=x/r and sin(theta)=y/r:
//  
//                  | x/r -y  0 |
//              J = | y/r  x  0 |
//                  | 0    0  1 |


// Initialize ptrList of cartesian fields
PtrList<volSymmTensorField> cartesianFields;

forAll(cylFields, i)
{
    // Get volSymmTensorField name
    word vstfName(cylFields[i].name());
    vstfName.replace("Cyl","");

    // Create volSymmTensorField
    volSymmTensorField vstf
    (
        IOobject
        (
            vstfName,
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("", cylFields[i].dimensions(), symmTensor::zero),
        "zeroGradient"
    );

    // Initialize Jacobian matrix
    tensor J(tensor::zero);
    
    forAll(mesh.C(), cellI)
    {
        // Get cell center
        vector Ci(mesh.C()[cellI]);

        // Compute radius
        const scalar r = Foam::sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1]);

        // Build Jacobian
        J.xx() = Ci[0]/r;
        J.xy() = -Ci[1];
        J.yx() = Ci[1]/r;
        J.yy() = Ci[0]; 
        J.zz() = 1;

        // Transform the field
        vstf[cellI] = symm ( J & cylFields[i][cellI] & J.T() ); // Check if use of symm is correct, not sure.
    }

    vstf.correctBoundaryConditions();

    // Write to disk
    vstf.write();

    Info << tab << cylFields[i].name() << " transformed into " << vstfName << "."
            << endl;

}