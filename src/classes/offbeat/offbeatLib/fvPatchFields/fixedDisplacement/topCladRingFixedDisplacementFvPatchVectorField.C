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

#include "topCladRingFixedDisplacementFvPatchVectorField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "globalOptions.H"
#include "gapGasModel.H"

namespace Foam
{
 
//* * * * * * * * * * * * * Protected member functions  * * * * * * * * * * *//


//  * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //
 
topCladRingFixedDisplacementFvPatchVectorField::topCladRingFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    coolantPressureFixedDisplacementFvPatchVectorField(p, iF),
    topCapRadius_(0.0),
    topCapPlenumRadius_(0.0),
    fuelTopPatchID_(-1),
    planeStrain_(true),
    displacementName_("D")
{}

topCladRingFixedDisplacementFvPatchVectorField::topCladRingFixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    coolantPressureFixedDisplacementFvPatchVectorField(p, iF, dict, false),
    topCapRadius_(readScalar(dict.lookup("topCapOuterRadius"))),
    topCapPlenumRadius_(readScalar(dict.lookup("topCapInnerRadius"))),
    fuelTopPatchID_(-1),
    planeStrain_(dict.lookupOrDefault<bool>("planeStrain", true)),
    displacementName_(dict.lookupOrDefault<word>("displacementName", "D"))
{}
 
 
topCladRingFixedDisplacementFvPatchVectorField::topCladRingFixedDisplacementFvPatchVectorField
(
    const topCladRingFixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    coolantPressureFixedDisplacementFvPatchVectorField(ptf, p, iF, mapper, mappingRequired),
    topCapRadius_(ptf.topCapRadius_),
    topCapPlenumRadius_(ptf.topCapPlenumRadius_),
    fuelTopPatchID_(ptf.fuelTopPatchID_),
    planeStrain_(true),
    displacementName_("D")
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
void Foam::topCladRingFixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Set fuelTopPatchID_
    if(fuelTopPatchID_ == -1)
    {
        const volVectorField& vfD =
            this->db().lookupObject<volVectorField>(displacementName_);

        int nFuelTopPatch(0);

        forAll(vfD.boundaryField(), patchID)
        {
            if(isA<plenumSpringBase>(vfD.boundaryField()[patchID]))
            {
                fuelTopPatchID_ = patchID;
                nFuelTopPatch += 1;
            }
        }

        if(nFuelTopPatch > 1)
        {
            FatalErrorInFunction()
            << "More than one plenumSpring type patch was found." << nl
            << "topCladRingPressure patchField requires one and "
            << "only one plenum spring patch-type in the model." 
            << abort(FatalError);
        }

        if(nFuelTopPatch == 0)
        {
            FatalErrorInFunction()
            << "No plenumSpring type patch was found." << nl
            << "topCladRingPressure patchField requires one and "
            << "only one plenum spring patch-type in the model." 
            << abort(FatalError);
        }
    }

    // Read angularFraction
    const globalOptions& globalOpt
    (patch().boundaryMesh().mesh().lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();

    // Include the gap gas pressure (only when the cap is not explicitly modeled)
    const gapGasModel& gapGas
    = db().lookupObject<gapGasModel>("gapGas");

    // Pi
    scalar pi = Foam::constant::mathematical::pi;

    // Update fluid pressure (if time list present)
    if(coolantPressureList_.valid())
    {
        coolantPressure_ = 
        coolantPressureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

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

    // Reference to total displacement field
    const volVectorField& D =
        db().lookupObject<volVectorField>(displacementName_);

    // Request spring elongation and spring modulus from fuel top patch
    const plenumSpringBase& fuelTopD(
        refCast<const plenumSpringBase>(D.boundaryField()[fuelTopPatchID_]));

    scalar Dtot(fuelTopD.springElongation());
    scalar springModulus(fuelTopD.springModulus());

    // Current average patch displacement
    const fvPatchField<vector>& Dp = D.boundaryField()[patch().index()]; 

    // Current internal displacement
    const scalarField Dz(Dp.patchInternalField() & n);

    // Top cap area fluid side. 
    scalar topCapArea = pi*pow(topCapRadius_, 2.0);

    // Top cap area plenum side. 
    scalar topCapPlenumArea = pi*pow(topCapPlenumRadius_, 2.0);

    // Explicit component of stress (sigma = twoMuLambda*grad + Q)
    const scalarField Q( ( (stress&n) - twoMuLambda*(gradField & n)) & n);

    // From the force balance:
    //         int(sigma_z dA) = k*D_s + p_g*Ai - p_c*Ao
    // where:
    //   - k is the spring modulus
    //   - D_s is the spring compression
    //   - p_g is the gas pressure acting on the inner cap area Ai
    //   - p_c is the coolant pressure acting on the outer cap area Ao

    // Prepare fields
    scalarField twoMuLambdaA(patch().size(), 0.0);
    scalarField twoMuLambdaDfA(patch().size(), 0.0);
    scalarField QA(patch().size(), 0.0);

    if(!planeStrain_)
    {
        // Assume plane stress --> local quantities:
        //         sigma_z * A 
        //      = (2mu+lambda)*A*gradD_z + Q*A = k*D_s + p_g*Ai - p_c*Ao
        //      = (2mu+lambda)*A*(Dfz - Dpz)/delta + Q*A = k*D_s + p_g*Ai - p_c*Ao
        //
        // Thus:
        //         Dfz = 1/((2mu+lambda)*A/delta)*[
        //         k*D_s + p_g*Ai - p_c*Ao + (2mu+lambda)*Dpz*A*/delta - Q*A]

        twoMuLambdaA = (twoMuLambda/delta)*gSum(magSf)/angularFraction;

        twoMuLambdaDfA = (twoMuLambda*Dz/delta)*gSum(magSf)/angularFraction;

        QA = Q*gSum(magSf)/angularFraction;
    }
    else
    {
        // Apply plane strain --> integral quantities:
        //         sigma_z * A 
        //      = int((2mu+lambda)*gradD_z*dA) + int(Q*dA) = k*D_s + p_g*Ai - p_c*Ao
        //      = int((2mu+lambda)*(Dfz - Dpz)/delta*dA) + int(Q*dA) = k*D_s + p_g*Ai - p_c*Ao
        //
        // Thus:
        //         Dfz = 1/int((2mu+lambda)*dA/delta)*[
        //         k*D_s + p_g*Ai - p_c*Ao + int((2mu+lambda)*Dpz*dA*/delta) - int(Q*dA)]

        // Calculate integrals
        twoMuLambdaA = gSum(twoMuLambda*magSf/delta/angularFraction);

        twoMuLambdaDfA = gSum(twoMuLambda*Dz*magSf/delta/angularFraction);

        QA = (gSum(Q*magSf/angularFraction));
    }     

    // Calculate spring and fluid force on the cladding
    scalar springForce = (springModulus*Dtot);
    scalarField fluidForce = (coolantPressure_*topCapArea);

    // Calculate plenum gas force
    // (positive because the plenum spring acts in the opposite direction on 
    // the cladding cap) 
    scalar plenumGasForce = gapGas.p()*topCapPlenumArea;

    // Set displacement
    vectorField disp = n*(springForce - fluidForce + plenumGasForce
        - QA + twoMuLambdaDfA)/twoMuLambdaA;

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

 
void Foam::topCladRingFixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    coolantPressureFixedDisplacementFvPatchVectorField::write(os);
#ifdef OPENFOAMFOUNDATION    

    writeEntry<scalar>(os, "topCapOuterRadius", topCapRadius_);
    writeEntry<scalar>(os, "topCapInnerRadius", topCapPlenumRadius_);

    writeEntryIfDifferent<bool>(os, "planeStrain", true, planeStrain_);

    writeEntryIfDifferent<word>(os, "displacementName", "D", displacementName_);
#elif OPENFOAMESI

    os.writeEntry<scalar>("topCapOuterRadius", topCapRadius_);
    os.writeEntry<scalar>("topCapInnerRadius", topCapPlenumRadius_);

    os.writeEntryIfDifferent<bool>("planeStrain", true, planeStrain_);

    os.writeEntryIfDifferent<word>("displacementName", "D", displacementName_);
#endif    
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    topCladRingFixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam