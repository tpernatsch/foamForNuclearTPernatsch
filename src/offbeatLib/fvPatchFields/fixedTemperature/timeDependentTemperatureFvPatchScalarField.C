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

#include "timeDependentTemperatureFvPatchScalarField.H"
#include "scalarFieldFieldINew.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

namespace Foam
{
 
//  * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
timeDependentTemperatureFvPatchScalarField::
timeDependentTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(p, iF),
    tempSeries_()
{}

timeDependentTemperatureFvPatchScalarField::
timeDependentTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fixedTemperatureFvPatchScalarField(p, iF, dict, valueRequired),    
    tempSeries_()
{          
#ifdef OPENFOAMFOUNDATION 
    tempSeries_.set
    ( 
        new Function1s::Table<scalar>
        (
            "temperatureSeries", dict.subDict("temperatureSeries")
        )
    );
#elif OPENFOAMESI
    tempSeries_.reset
    ( 
        new Table
        (
            "temperatureSeries", dict.subDict("temperatureSeries")
        )
    );
#endif

    T_ = 
    tempSeries_->value
    (
        patch().boundaryMesh().mesh().time().timeToUserTime
        (
            this->db().time().value()
        )
    );
 }
 
 
timeDependentTemperatureFvPatchScalarField::
timeDependentTemperatureFvPatchScalarField
(
    const timeDependentTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedTemperatureFvPatchScalarField(ptf, p, iF, mapper),
    tempSeries_()
{}
 
 
timeDependentTemperatureFvPatchScalarField::
timeDependentTemperatureFvPatchScalarField
(
    const timeDependentTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedTemperatureFvPatchScalarField(ptf, iF),
    tempSeries_()
{}
 
 
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 

void Foam::timeDependentTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    T_ = 
    tempSeries_->value
    (
        patch().boundaryMesh().mesh().time().timeToUserTime
        (
            this->db().time().value()
        )
    );

    fixedTemperatureFvPatchScalarField::updateSurfaceTemperature();

    fvPatchField<scalar>::updateCoeffs();
}
 
void Foam::timeDependentTemperatureFvPatchScalarField::write(Ostream& os) const
{
    fixedTemperatureFvPatchScalarField::write(os);

    if (tempSeries_.valid())
    {
        os.writeKeyword("temperatureSeries") << nl;
        os << token::BEGIN_BLOCK << nl;
#ifdef OPENFOAMFOUNDATION  
        tempSeries_->write(os);
#elif OPENFOAMESI
        tempSeries_->writeData(os);
#endif
        os << token::END_BLOCK << nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    timeDependentTemperatureFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam