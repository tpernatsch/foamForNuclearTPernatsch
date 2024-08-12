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

#include "YoungModulusMolybdenum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusMolybdenum, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel, 
        YoungModulusMolybdenum, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusMolybdenum::YoungModulusMolybdenum
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, defaultModel),       
    par1(3.349E+11),
    par2(5.101E+07),
    perturb(1.0)
{
    if(dict.found("YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.subDict("YoungModulus");

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 3.349E+11);
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 5.101E+07);

        perturb = YoungModulusDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusMolybdenum::~YoungModulusMolybdenum()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusMolybdenum::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{    

    forAll(addr, i)
    {
        const label cellI = addr[i];
        const scalar nominalValue = 3.349E+11 - 5.101E+07*T[cellI];

        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
