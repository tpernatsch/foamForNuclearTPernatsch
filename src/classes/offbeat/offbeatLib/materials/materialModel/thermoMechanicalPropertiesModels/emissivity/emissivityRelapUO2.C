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

#include "emissivityRelapUO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(emissivityRelapUO2, 0);
    addToRunTimeSelectionTable
    (
        emissivityModel, 
        emissivityRelapUO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emissivityRelapUO2::emissivityRelapUO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    emissivityModel(mesh, dict, defaultModel),
    par1(0.7856),
    par2(1.5263e-5),
    perturb(1.0)
{
    if(dict.found("emissivity"))
    {
        const dictionary& emissivityDict = dict.subDict("emissivity");

        par1 = emissivityDict.lookupOrDefault<scalar>("par1", 0.7856);
        par2 = emissivityDict.lookupOrDefault<scalar>("par2", 1.5263e-5);

        perturb = emissivityDict.lookupOrDefault<scalar>("perturb", 1.0);
        
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::emissivityRelapUO2::~emissivityRelapUO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::emissivityRelapUO2::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{  
    forAll(addr, i)
    {   
        const label cellI = addr[i];      

        const scalar Ti = T[cellI];
        
        const scalar nominalValue =  par1 + par2*Ti;

        sf[cellI] = nominalValue*perturb;
    }
}

    
 
// ************************************************************************* //
