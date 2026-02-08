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

#include "fixedDisplacementZeroShearFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "transformField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_()
{
    w_ = 0;
}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    solidDirectionMixedFvPatchVectorField(ptf, p, iF, mapper),
    totalDisp_(p.size()),
    dispSeries_()
{
    w_ = 0;

    if (mappingRequired)
    {
        if (mapper.hasUnmapped())
        {
         totalDisp_ = vector::zero;
        }

    #ifdef OPENFOAMFOUNDATION
        mapper(totalDisp_, ptf.totalDisp_);
    #elif OPENFOAMESI
        totalDisp_.autoMap(mapper);
    #endif
    }
}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    solidDirectionMixedFvPatchVectorField(p, iF),
    totalDisp_(p.size(), Foam::vector(0, 0, 0)),
    dispSeries_()
{
    if(valueRequired)
    {
        // Check if displacement is time-varying
        if (dict.found("displacementSeries"))
        {
#ifdef OPENFOAMFOUNDATION            
            dispSeries_.set
            ( 
                new Function1s::Table<vector>
                (
                    "displacementSeries", dict.subDict("displacementSeries")
                )
            );
#elif OPENFOAMESI            
            dispSeries_.reset
            ( 
                Function1<vector>::New
                (
                    "displacementSeries", dict.subDict("displacementSeries")
                )
            );
#endif  

            refValue() = 
            dispSeries_->value
            (
                patch().boundaryMesh().mesh().time().timeToUserTime
                (
                    this->db().time().value()
                )
            );
        }
        else if (dict.found("value"))
        {
            totalDisp_ = vectorField("value", dict, p.size());
            refValue() = totalDisp_;
        }
        else
        {
            FatalErrorIn
            (
                "fixedDisplacementZeroShearFvPatchVectorField::"
                "fixedDisplacementZeroShearFvPatchVectorField"
            )   << "displacementSeries nor value entries not found for patch " 
                << patch().name()
                << abort(FatalError);
        }    
    }

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
    
    // Set w_ to zero (for compatibility with all other directionMixed BC)
    w_ = 0; 

    this->valueFraction() = sqr(patch().nf());

    Field<vector> normalValue = transform(valueFraction(), refValue());

    Field<vector> gradValue =
        this->patchInternalField() + refGrad()/this->patch().deltaCoeffs();

    const symmTensorField valueFractionCopy = valueFraction();
    Field<vector> transformGradValue = 
        transform(I - valueFractionCopy, gradValue);

    Field<vector>::operator=(normalValue + transformGradValue);
}


fixedDisplacementZeroShearFvPatchVectorField::
fixedDisplacementZeroShearFvPatchVectorField
(
    const fixedDisplacementZeroShearFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    solidDirectionMixedFvPatchVectorField(ptf, iF),
    totalDisp_(),
    dispSeries_()
{
    w_ = 0;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void fixedDisplacementZeroShearFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    solidDirectionMixedFvPatchVectorField::autoMap(m);
#ifdef OPENFOAMFOUNDATION
    m(totalDisp_, totalDisp_);
#elif OPENFOAMESI
    totalDisp_.autoMap(m);
#endif
}


// Reverse-map the given fvPatchField onto this fvPatchField
void fixedDisplacementZeroShearFvPatchVectorField::rmap
(
    const fvPatchField<vector>& ptf,
    const labelList& addr
)
{
    solidDirectionMixedFvPatchVectorField::rmap(ptf, addr);

    const fixedDisplacementZeroShearFvPatchVectorField& dmptf =
        refCast<const fixedDisplacementZeroShearFvPatchVectorField>(ptf);

    totalDisp_.rmap(dmptf.totalDisp_, addr);
}


void fixedDisplacementZeroShearFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    // Find whether or not the simulation usese a small-strain incremental
    // solver with mesh update, to apply a special treatement
    bool smallStrainIncUpdated = false;
    if( (internalField().name() == "DD") and 
        not(patch().boundaryMesh().mesh().foundObject<volTensorField>("Finv")))
    {
        smallStrainIncUpdated = true;
    }

    vectorField disp = totalDisp_;

    if (dispSeries_.valid())
    {
        disp = 
        dispSeries_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    if (internalField().name() == "DD")
    {
        // Incremental approach, so we wil set the increment of displacement
        // Lookup the old displacement field and subtract it from the total
        // displacement
        const volVectorField& Dold =
            db().lookupObject<volVectorField>("D").oldTime();

        disp -= Dold.boundaryField()[patch().index()];
    }

    // Set displacement
    refValue() = disp;

    // Set gradient to force zero shear traction
    const fvPatchField<scalar>& mu =
        patch().lookupPatchField<volScalarField, scalar>("mu");

    const fvPatchField<scalar>& lambda =
        patch().lookupPatchField<volScalarField, scalar>("lambda");
        
    const fvPatchField<tensor>& sigmaExp = smallStrainIncUpdated ?
        patch().lookupPatchField<volTensorField, tensor>("DSigmaExp") :
        patch().lookupPatchField<volTensorField, tensor>("sigmaExp");
    
    scalarField twoMuLambda(2*mu + lambda);
    
    vectorField n(patch().nf());  
    
    refGrad() = -(n & sigmaExp) / twoMuLambda ;

    solidDirectionMixedFvPatchVectorField::updateCoeffs();
}


// Write
void fixedDisplacementZeroShearFvPatchVectorField::write(Ostream& os) const
{
    solidDirectionMixedFvPatchVectorField::write(os);

    if (dispSeries_.valid())
    {
        os.writeKeyword("displacementSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
#ifdef OPENFOAMFOUNDATION        
        dispSeries_->write(os);
#elif OPENFOAMESI        
        dispSeries_->writeData(os);
#endif        
        os << token::END_BLOCK << nl;
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedDisplacementZeroShearFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
