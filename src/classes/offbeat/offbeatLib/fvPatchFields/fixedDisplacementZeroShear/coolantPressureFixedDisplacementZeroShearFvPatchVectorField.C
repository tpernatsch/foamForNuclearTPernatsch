/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Hrvoje Jasak
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "coolantPressureFixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coolantPressureFixedDisplacementZeroShearFvPatchVectorField::
coolantPressureFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementZeroShearFvPatchVectorField(p, iF),
    coolantPressure_(p.size(), 0.0),
    coolantPressureList_(),
    planeStrain_(true),
    displacementName_("D")
{}


coolantPressureFixedDisplacementZeroShearFvPatchVectorField::
coolantPressureFixedDisplacementZeroShearFvPatchVectorField
(
    const coolantPressureFixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedDisplacementZeroShearFvPatchVectorField(ptf, p, iF, mapper),
    coolantPressure_(ptf.coolantPressure_),
    coolantPressureList_(),
    planeStrain_(true),
    displacementName_("D")
{}


coolantPressureFixedDisplacementZeroShearFvPatchVectorField::
coolantPressureFixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedDisplacementZeroShearFvPatchVectorField(p, iF, dict, false),
    coolantPressure_(p.size(), 0),
    coolantPressureList_(),
    planeStrain_(dict.lookupOrDefault<bool>("planeStrain", true)),
    displacementName_(dict.lookupOrDefault<word>("displacementName", "D"))
{
    //- Read fluid pressure (either list or fixed pressure)
    if(dict.found("coolantPressureList"))
    {
#ifdef OPENFOAMFOUNDATION            
        coolantPressureList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "coolantPressureList", dict.subDict("coolantPressureList")
            )
        );
#elif OPENFOAMESI            
        coolantPressureList_.reset
        ( 
            Function1<scalar>::New
            (
                "coolantPressureList", dict.subDict("coolantPressureList")
            )
        );
#endif  
        
        coolantPressure_ = 
        coolantPressureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }
    else if (dict.found("coolantPressure"))
    {
        coolantPressure_ = scalarField("coolantPressure", dict, p.size());
    }
    else
    {
        FatalErrorInFunction() 
        << "elasticSpring BC for " << patch().name() << " requires:" << nl
        << "- \"coolantPressureList\" or" << nl
        << "- \"coolantPressure\"" << abort(FatalError) << endl;
    }

    refValue() = vectorField("value", dict, p.size());
    
    //- Read gradient if present (otherwise zero gradient)
    if (dict.found("refGradient"))
    {
        refGrad() =
        (
            vectorField("refGradient", dict, p.size())
        );
    }
    else
    {
        refGrad() = vector::zero;
    }  

    this->valueFraction() = sqr(patch().nf());

    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    Field<vector> transformGradValue =
        transform(I - valueFraction(), gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);
}


coolantPressureFixedDisplacementZeroShearFvPatchVectorField::
coolantPressureFixedDisplacementZeroShearFvPatchVectorField
(
    const coolantPressureFixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedDisplacementZeroShearFvPatchVectorField(ptf, iF),
    coolantPressure_(ptf.coolantPressure_),
    coolantPressureList_(),
    planeStrain_(true),
    displacementName_("D")
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void coolantPressureFixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const globalOptions& globalOpt
    (patch().boundaryMesh().mesh().lookupObject<globalOptions>("globalOptions"));

    scalar angularFraction = globalOpt.angularFraction();
    vector pinDirection = globalOpt.pinDirection();

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

    // Average patch normal
    vector nAvg = gSum(n*magSf)/gSum(magSf);

    // Current average patch displacement
    const fvPatchField<vector>& D =
        patch().lookupPatchField<volVectorField, vector>(displacementName_);   

    // Current internal displacement
    const scalarField Dz(D.patchInternalField() & n);

    // Explicit component of stress (sigma = twoMuLambda*grad + Q)
    const scalarField Q( ( (stress&n) - twoMuLambda*(gradField & n)) & n);

    // From the force balance:
    //         int(sigma_z dA) = +-p_c*A
    // where:
    //   - k is the spring modulus
    //   - D_s is the spring compression
    //   - p_c is the coolant pressure
    //   - A is the area of this patch divided by the angular fraction

    // Prepare fields
    scalarField twoMuLambdaA(patch().size(), 0.0);
    scalarField twoMuLambdaDfA(patch().size(), 0.0);
    scalarField QA(patch().size(), 0.0);

    if(!planeStrain_)
    {
        // Assume plane stress --> local quantities:
        //         sigma_z * A 
        //      = (2mu+lambda)*A*gradD_z + Q*A = +-p_c*A
        //      = (2mu+lambda)*A*(Dfz - Dpz)/delta + Q*A = +-p_c*A
        //
        // Thus:
        //         Dfz = 1/((2mu+lambda)*A/delta)*[
        //         +-p_c*A + (2mu+lambda)*Dpz*A*/delta - Q*A]
        //
        // (sign of gradient and of forces change depending on the patch being 
        // top cap inner or bottom cap inner)

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
        //      = int((2mu+lambda)*gradD_z*dA) + int(Q*dA) = +-p_c*A
        //      = int((2mu+lambda)*(Dfz - Dpz)/delta*dA) + int(Q*dA) = +-p_c*A
        //
        // Thus:
        //         Dfz = 1/int((2mu+lambda)*dA/delta)*[
        //         +-p_c*A + int((2mu+lambda)*Dpz*dA*/delta) - int(Q*dA)]
        //
        // (sign of gradient and of forces change depending on the patch being 
        // top cap inner or bottom cap inner)

        // Calculate integrals
        twoMuLambdaA = (sign(nAvg & pinDirection)*(
            gSum(twoMuLambda*magSf/delta/angularFraction)));

        twoMuLambdaDfA = (sign(nAvg & pinDirection)*(
            gSum(twoMuLambda*Dz*magSf/delta/angularFraction)));

        QA = (gSum(Q*magSf/angularFraction));
    }     

    // Calculate fluid force
    scalar fluidForce = 
    sign(nAvg & pinDirection)*gSum(coolantPressure_*magSf)/angularFraction;

    // Set displacement
    vectorField disp
    (
        pinDirection*(-fluidForce - QA + twoMuLambdaDfA)/twoMuLambdaA
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

    refValue() = disp;  
    
    // Set gradient
    refGrad() = 
    (
        - (n & (stress - gradField*twoMuLambda))
    ) / twoMuLambda ;

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void coolantPressureFixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    fixedDisplacementZeroShearFvPatchVectorField::write(os);
    
#ifdef OPENFOAMFOUNDATION
    if (coolantPressureList_.valid())
    {
        os.writeKeyword("coolantPressureList") << nl;
        os << token::BEGIN_BLOCK << nl;
        coolantPressureList_->write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        writeEntry<scalarField>(os, "coolantPressure", coolantPressure_);
    }

    writeEntryIfDifferent<bool>(os, "planeStrain", true, planeStrain_);
    writeEntryIfDifferent<word>(os, "displacementName", "D", displacementName_);
#elif OPENFOAMESI
    if (coolantPressureList_.valid())
    {
        os.writeKeyword("coolantPressureList") << nl;
        os << token::BEGIN_BLOCK << nl;
        coolantPressureList_->writeData(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        coolantPressure_.writeEntry( "coolantPressure", os);
    }

    os.writeEntryIfDifferent<bool>("planeStrain", true, planeStrain_);

    os.writeEntryIfDifferent<word>("displacementName", "D", displacementName_);
#endif      
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    coolantPressureFixedDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
