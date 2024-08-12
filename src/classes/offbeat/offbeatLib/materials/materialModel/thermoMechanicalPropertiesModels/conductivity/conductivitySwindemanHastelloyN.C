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

#include "conductivitySwindemanHastelloyN.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivitySwindemanHastelloyN, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel,
        conductivitySwindemanHastelloyN,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivitySwindemanHastelloyN::conductivitySwindemanHastelloyN
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, defaultModel),
    par1(8.431),
    par2(0.0205),
    perturb(1)
{
    if(dict.found("conductivity"))
    {
        const dictionary& conductivityDict = dict.subDict("conductivity");

        par1 = conductivityDict.lookupOrDefault<scalar>("par1", 8.431);
        par2 = conductivityDict.lookupOrDefault<scalar>("par2", 0.0205);

        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivitySwindemanHastelloyN::~conductivitySwindemanHastelloyN()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivitySwindemanHastelloyN::correct
(
    scalarField& sf,
    const scalarField& T,
    const labelList& addr
)
{

    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T[cellI];

        const scalar nominalValue =  par1 + par2*(Ti-273.15) ;

        sf[cellI] = nominalValue*perturb;
    }
}


    


// ************************************************************************* //

