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

#include "fixedDisplacementFvPatchVectorField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{
 
//* * * * * * * * * * * * * Protected member functions  * * * * * * * * * * *//


//  * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //
 
fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    totalDisp_(p.size(), vector::zero),
    dispSeries_()
{}

fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    totalDisp_(p.size(), vector::zero),
    dispSeries_()
{
    if(valueRequired)
    {
        // Check if displacement is time-varying
        if (dict.found("displacementSeries") and internalField().name() != "DD")
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

            totalDisp_ = 
            dispSeries_->value
            (
                patch().boundaryMesh().mesh().time().timeToUserTime
                (
                    this->db().time().value()
                )
            );
            
            fvPatchField<vector>::operator==(totalDisp_);
        }
        else if (dict.found("value"))
        {
            totalDisp_ =  vectorField("value", dict, p.size());
            
            fvPatchField<vector>::operator==(totalDisp_);
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

}
 
 
fixedDisplacementFvPatchVectorField::fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
#ifdef OPENFOAMFOUNDATION
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper, mappingRequired),
#elif OPENFOAMESI
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
#endif        
    totalDisp_(ptf.totalDisp_),
    dispSeries_()
{}
 
 
fixedDisplacementFvPatchVectorField::
fixedDisplacementFvPatchVectorField
(
    const fixedDisplacementFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    totalDisp_(ptf.totalDisp_),
    dispSeries_()
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
void Foam::fixedDisplacementFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vectorField currentDisp = totalDisp_;

    if (dispSeries_.valid())
    {
        currentDisp = 
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

        currentDisp -= Dold.boundaryField()[patch().index()];
    }

    fvPatchField<vector>::operator==(currentDisp);

    fvPatchField<vector>::updateCoeffs();
}

 
void Foam::fixedDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchField<vector>::write(os);

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
    fixedDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam