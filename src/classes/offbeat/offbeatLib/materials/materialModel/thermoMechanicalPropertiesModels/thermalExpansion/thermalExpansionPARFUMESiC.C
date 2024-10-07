/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "thermalExpansionPARFUMESiC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionPARFUMESiC, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel,
        thermalExpansionPARFUMESiC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionPARFUMESiC::thermalExpansionPARFUMESiC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    alpha_(dict.lookupOrDefault("thermalExpansionCoefficient", 4.9e-6))
{
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionPARFUMESiC::~thermalExpansionPARFUMESiC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionPARFUMESiC::correct
(
    symmTensorField& sf,
    const scalarField& T,
    const labelList& addr
)
const
{
    forAll(addr, i)
    {
        const label cellI = addr[i];

        // Temperature
        const scalar Ti = T[cellI];
        const scalar Tref = Tref_.value();

        symmTensor nominalValue = alpha_ * (Ti-Tref)*I;


        // NOTE: avoids instabilities when trying to simulate a material at
        // constant temperature equal to Tref
        if
        (
            nominalValue.xx() < 1e-7 &&
            nominalValue.yy() < 1e-7 &&
            nominalValue.zz() < 1e-7
        )
        {
            nominalValue *= 0;
        }

        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
