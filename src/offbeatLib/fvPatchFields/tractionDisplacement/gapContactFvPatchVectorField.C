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

#include "gapContactFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//
    

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    contactFvPatchVectorField(p, iF)
{}


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const gapContactFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    contactFvPatchVectorField(ptf, p, iF, mapper)
{}


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    contactFvPatchVectorField(p, iF, dict)
{}


gapContactFvPatchVectorField::gapContactFvPatchVectorField
(
    const gapContactFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    contactFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



void gapContactFvPatchVectorField::updateTraction()
{
    contactFvPatchVectorField::updateTraction();
    
    const objectRegistry& db = this->patch().boundaryMesh().mesh();

    // Include the gap gas pressure
    const gapGasModel& gapGas = db.lookupObject<gapGasModel>("gapGas");

    pressure_ += gapGas.p();
    
    // TODO: 
    // - should the gapGas pressure disappear?
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    gapContactFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
