/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "gapPressureFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

gapPressureFvPatchVectorField::
gapPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF)
{}


gapPressureFvPatchVectorField::
gapPressureFvPatchVectorField
(
    const gapPressureFvPatchVectorField& fplpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(fplpvf, p, iF, mapper)
{}


gapPressureFvPatchVectorField::
gapPressureFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict, false)
{
    // Read initial pressure value if present (for accurate restart)
    if(dict.found("pressure"))
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }

}


gapPressureFvPatchVectorField::
gapPressureFvPatchVectorField
(
    const gapPressureFvPatchVectorField& fplpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(fplpvf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void gapPressureFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


void gapPressureFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void gapPressureFvPatchVectorField::updateTraction()
{    
    // Include the gap gas pressure
    const gapGasModel& gapGas
        = patch().boundaryMesh().mesh().lookupObject<gapGasModel>("gapGas");

    pressure_ = gapGas.p();
    traction_ = vector::zero;    
}


void gapPressureFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    gapPressureFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
