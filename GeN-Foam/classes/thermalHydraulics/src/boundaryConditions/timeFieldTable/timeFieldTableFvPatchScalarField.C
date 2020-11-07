/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "timeFieldTableFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    currentIndex_(0),
    table_(),
    tStart_(0)
{}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict, false),
    currentIndex_(0),
    table_
    (
        dict.get
        <
            List
            <
                Tuple2
                <
                    scalar, 
                    scalarField
                >
            >
        >("table")
    ),
    tStart_(table_[0].first())
{    
    updateField();
}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const timeFieldTableFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    currentIndex_(ptf.currentIndex_),
    table_(ptf.table_),
    tStart_(ptf.tStart_)
{}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const timeFieldTableFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    currentIndex_(ptf.currentIndex_),
    table_(ptf.table_),
    tStart_(ptf.tStart_)
{}


Foam::timeFieldTableFvPatchScalarField::
timeFieldTableFvPatchScalarField
(
    const timeFieldTableFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    currentIndex_(ptf.currentIndex_),
    table_(ptf.table_),
    tStart_(ptf.tStart_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeFieldTableFvPatchScalarField::updateField()
{
    //- Simulation time
    const scalar t = db().time().timeOutputValue();

    if (t > tStart_)
    {
        scalar t0(0);
        scalar t1(0);
        if (currentIndex_ < (table_.size()-1))
        {
            t0 = table_[currentIndex_].first();
            t1 = table_[currentIndex_+1].first();
            while (true)
            {
                if ((t >= t0 and t <= t1))
                    break;
                currentIndex_++;
                t0 = table_[currentIndex_].first();
                t1 = table_[currentIndex_+1].first(); 
            }
        }

        if (currentIndex_ < (table_.size()-1))
        {
            //- Read table values at correct index
            scalarField& f0(table_[currentIndex_].second());
            scalarField& f1(table_[currentIndex_+1].second());

            //- Interpolation coefficient
            scalar c((t1-t)/(t1-t0));

            //- Linear interpolation between the fields at the provided times
            //  if the time falls in between two time bins 
            this->operator==
            (
                c*f0+(1.0-c)*f1
            );
        }
        else //- i.e. t > tLast
        {
            this->operator==(table_[table_.size()-1].second());
        }
    }
    else //- i.e. t <= tStart
    {
        this->operator==(table_[0].second());
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::timeFieldTableFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateField();
}


void Foam::timeFieldTableFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry
    <
        List
        <
            Tuple2
            <
                scalar, 
                scalarField
            >
        >
    >("table", table_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchScalarField,
       timeFieldTableFvPatchScalarField
   );
}


// ************************************************************************* //
