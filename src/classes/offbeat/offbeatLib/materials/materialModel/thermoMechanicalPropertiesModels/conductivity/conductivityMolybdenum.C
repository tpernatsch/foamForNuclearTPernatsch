/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022 OpenFOAM Foundation
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

#include "conductivityMolybdenum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityMolybdenum, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityMolybdenum, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityMolybdenum::conductivityMolybdenum
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    conductivityModel(mesh, dict, defaultModel),       
    par1(9.128E-06),
    par2(4.945E-2),
    par3(152),
    perturb(1)
{
    if(dict.found("conductivity"))
    {
        const dictionary& conductivityDict = dict.subDict("conductivity");

        par1 = conductivityDict.lookupOrDefault<scalar>("par1", 9.128E-06);
        par2 = conductivityDict.lookupOrDefault<scalar>("par2", 4.945E-2);
        par3 = conductivityDict.lookupOrDefault<scalar>("par3", 152);

        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityMolybdenum::~conductivityMolybdenum()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityMolybdenum::correct
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
        
        const scalar nominalValue = par1*Ti*Ti - par2*Ti + par3;

        sf[cellI] = nominalValue*perturb;
    }
}

// ************************************************************************* //
