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

#include "velocityRundownFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "one.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityRundownFvPatchVectorField::
velocityRundownFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U0_(p.size(), vector::zero),
    A_(),
    B_(),
    C_(),
    t0_()
{}


Foam::velocityRundownFvPatchVectorField::
velocityRundownFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    U0_("value", dict, p.size()),
    A_(dict.get<scalar>("const")),
    B_(dict.get<scalar>("coeff")),
    C_(dict.get<scalar>("exp")),
    t0_(dict.get<scalar>("startTime"))
{
    const scalar t = db().time().timeOutputValue();
    if (t>t0_)
    {
        U0_ /= Foam::pow(A_+B_*(t-t0_), C_);
    }

    this->operator==(U0_);
}


Foam::velocityRundownFvPatchVectorField::
velocityRundownFvPatchVectorField
(
    const velocityRundownFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    U0_(ptf.U0_),
    A_(ptf.A_),
    B_(ptf.B_),
    C_(ptf.C_),
    t0_(ptf.t0_)
{}


Foam::velocityRundownFvPatchVectorField::
velocityRundownFvPatchVectorField
(
    const velocityRundownFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    U0_(ptf.U0_),
    A_(ptf.A_),
    B_(ptf.B_),
    C_(ptf.C_),
    t0_(ptf.t0_)
{}


Foam::velocityRundownFvPatchVectorField::
velocityRundownFvPatchVectorField
(
    const velocityRundownFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    U0_(ptf.U0_),
    A_(ptf.A_),
    B_(ptf.B_),
    C_(ptf.C_),
    t0_(ptf.t0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityRundownFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();

    if (t>t0_)
    {
        vectorField Up(this->patch().size(), vector::zero);

        Up = U0_*Foam::pow(A_+B_*(t-t0_), C_);
        
        this->operator==(Up);

        fixedValueFvPatchVectorField::updateCoeffs();
    }
}


void Foam::velocityRundownFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry<scalar>("const", A_);
    os.writeEntry<scalar>("coeff", B_);
    os.writeEntry<scalar>("exp", C_);
    os.writeEntry<scalar>("startTime", t0_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       velocityRundownFvPatchVectorField
   );
}


// ************************************************************************* //
