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

#include "thermalExpansionPARFUMEBuffer.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionPARFUMEBuffer, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel,
        thermalExpansionPARFUMEBuffer,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionPARFUMEBuffer::thermalExpansionPARFUMEBuffer
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    par1(5.0),
    par2(0.11),
    par3(400.0),
    par4(700.0),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        par1 = thermalExpansionDict.lookupOrDefault<scalar>("par1", 5.0);
        par2 = thermalExpansionDict.lookupOrDefault<scalar>("par2", 0.11);
        par3 = thermalExpansionDict.lookupOrDefault<scalar>("par3", 400.0);
        par4 = thermalExpansionDict.lookupOrDefault<scalar>("par4", 700.0);

        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionPARFUMEBuffer::~thermalExpansionPARFUMEBuffer()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionPARFUMEBuffer::correct
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

        // Temperature must be given in °C
        const scalar Ti = T[cellI] - 273.15;
        const scalar Tref = Tref_.value() -273.15;

        const scalar alpha = par1*(1.0+par2*(Ti-par3)/par4)*1e-6;

        symmTensor nominalValue = alpha*(Ti-Tref)*I;

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

        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
