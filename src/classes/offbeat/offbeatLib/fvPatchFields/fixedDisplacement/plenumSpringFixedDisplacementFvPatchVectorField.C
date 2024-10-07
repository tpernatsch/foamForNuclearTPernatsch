/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "plenumSpringFixedDisplacementFvPatchVectorField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "gapGasModel.H"
#include "globalOptions.H"

namespace Foam
{
 
//* * * * * * * * * * * * * Protected member functions  * * * * * * * * * * *//


//  * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //
 
plenumSpringFixedDisplacementFvPatchVectorField::plenumSpringFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementFvPatchVectorField(p, iF),
    plenumSpringBase(p),
    planeStrain_(true)
{}

plenumSpringFixedDisplacementFvPatchVectorField::plenumSpringFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedDisplacementFvPatchVectorField(p, iF, dict, false),
    plenumSpringBase(p, dict),
    planeStrain_(dict.lookupOrDefault<bool>("planeStrain", true))
{
    // Set initial value
    vectorField disp(vectorField("value", dict, p.size()));
    
    fvPatchField<vector>::operator==(disp);
}
 
 
plenumSpringFixedDisplacementFvPatchVectorField::plenumSpringFixedDisplacementFvPatchVectorField
(
    const plenumSpringFixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fixedDisplacementFvPatchVectorField(ptf, p, iF, mapper, mappingRequired),
    plenumSpringBase(p, *this),
    planeStrain_(true)
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
void Foam::plenumSpringFixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //- Include the gap gas pressure
    const gapGasModel& gapGas
    = patch().boundaryMesh().mesh().lookupObject<gapGasModel>("gapGas");

    const globalOptions& globalOpt
    (patch().boundaryMesh().mesh().lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();
    vector pinDirection = globalOpt.pinDirection();

    // Lookup material properties and set twoMuLambda
    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
    
    scalarField twoMuLambda(2*mu + lambda);

    // Lookup boundary stress       
    const fvPatchField<symmTensor>& stress =
        patch().lookupPatchField<volSymmTensorField, symmTensor>("sigma");

    // Lookup the gradient field
    const fvPatchField<tensor>& gradField =
        patch().lookupPatchField<volTensorField, tensor>
        (
            "grad" + internalField().name()
        );

    // Polypatch quantities
    vectorField n(patch().nf());      
    scalarField delta(patch().delta()&n);  
    scalarField magSf(patch().magSf());  

    // Average patch normal
    vector nAvg = gSum(n*magSf)/gSum(magSf);
        
    // Current average patch displacement
    const fvPatchField<vector>& Dp = 
        patch().lookupPatchField<volVectorField, vector>(displacementName_);

    // Current internal displacement
    const scalarField Dz(Dp.patchInternalField() & n);

    // It is assumed that:
    //         sigma_z = (2mu+lambda)*gradD_z + Q

    // Explicit component of stress or Q
    const scalarField Q( ( (stress&n) - twoMuLambda*(gradField & n)) & n);

    // From the force balance:
    //         int(sigma_z dA) = -k*D_s - p_g*A
    // where:
    //   - k is the spring modulus
    //   - D_s is the spring compression
    //   - p_g is the gas pressure acting
    //   - A is the area of this patch divided by the angular fraction

    // Prepare fields
    scalarField twoMuLambdaA(patch().size(), 0.0);
    scalarField twoMuLambdaDfA(patch().size(), 0.0);
    scalarField QA(patch().size(), 0.0);

    if(!planeStrain_)
    {
        // Assume plane stress --> local quantities:
        //         sigma_z * A 
        //      = (2mu+lambda)*A*gradD_z + Q*A = -k*D_s - p_g*A
        //      = (2mu+lambda)*A*(Dfz - Dpz)/delta + Q*A = -k*D_s - p_g*A
        //
        // Thus:
        //         Dfz = 1/((2mu+lambda)*A/delta)*[
        //         -k*D_s - p_g*A + (2mu+lambda)*Dpz*A*/delta - Q*A]
        //
        // (sign of gradient and of forces change depending on the patch being 
        // fuel top or top cap inner)

        twoMuLambdaA = sign(nAvg & pinDirection)*
            (twoMuLambda/delta)*gSum(magSf)/angularFraction;

        twoMuLambdaDfA = sign(nAvg & pinDirection)*
            (twoMuLambda*Dz/delta)*gSum(magSf)/angularFraction;

        QA = Q*gSum(magSf)/angularFraction;
    }
    else
    {
        // Apply plane strain --> integral quantities:
        //         sigma_z * A 
        //      = int((2mu+lambda)*gradD_z*dA) + int(Q*dA) = -k*D_s - p_g*A
        //      = int((2mu+lambda)*(Dfz - Dpz)/delta*dA) + int(Q*dA) = -k*D_s - p_g*A
        //
        // Thus:
        //         Dfz = 1/int((2mu+lambda)*dA/delta)*[
        //         -k*D_s - p_g*A + int((2mu+lambda)*Dpz*dA*/delta) - int(Q*dA)]
        //
        // (sign of gradient and of forces change depending on the patch being 
        // fuel top or top cap inner)

        // Calculate integrals
        twoMuLambdaA = (sign(nAvg & pinDirection)*(
            gSum(twoMuLambda*magSf/delta/angularFraction)));

        twoMuLambdaDfA = (sign(nAvg & pinDirection)*(
            gSum(twoMuLambda*Dz*magSf/delta/angularFraction)));

        QA = (gSum(Q*magSf/angularFraction));
    }

    // Calculate total spring elongation (positive if spring contracts)
    scalar Dtot(plenumSpringBase::springElongation());

    // Calculate spring and gas force (sign changes based on patch being fuel
    // or topCap)
    scalar springForce = sign(nAvg & pinDirection)*springModulus_*Dtot;

    scalar gasForce    = 
    sign(nAvg & pinDirection)*gapGas.p()*gSum(magSf)/angularFraction;

    // Set displacement 
    vectorField disp
    (   
        pinDirection*(-springForce - gasForce - QA + 
            twoMuLambdaDfA)/twoMuLambdaA
    );

    if (internalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>(displacementName_).oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    fvPatchField<vector>::operator==(disp);

    fvPatchField<vector>::updateCoeffs();
}

 
void Foam::plenumSpringFixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedDisplacementFvPatchVectorField::write(os);
    plenumSpringBase::write(os);
#ifdef OPENFOAMFOUNDATION
    writeEntryIfDifferent<bool>(os, "planeStrain", true, planeStrain_);
#elif OPENFOAMESI
    os.writeEntryIfDifferent<bool>("planeStrain", true, planeStrain_);
#endif
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    plenumSpringFixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam