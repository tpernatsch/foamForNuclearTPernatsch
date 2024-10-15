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

#include "heatCapacityMolybdenum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityMolybdenum, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityMolybdenum, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityMolybdenum::heatCapacityMolybdenum
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, defaultModel),
    par1(9.74E-06),
    par2(5.37E-2),
    par3(235),
    perturb(1.0)
{
    if(dict.found("heatCapacity"))
    {
        const dictionary& heatCapacityDict = dict.subDict("heatCapacity");

        par1 = heatCapacityDict.lookupOrDefault<scalar>("par1", 9.74E-06);
        par2 = heatCapacityDict.lookupOrDefault<scalar>("par2", 5.37E-2);
        par3 = heatCapacityDict.lookupOrDefault<scalar>("par3", 235);

        perturb = heatCapacityDict.lookupOrDefault<scalar>("perturb", 1.0);
        
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityMolybdenum::~heatCapacityMolybdenum()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacityMolybdenum::correct
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
        
        const scalar nominalValue = par1*Ti*Ti + par2*Ti + par3;

        sf[cellI] = nominalValue*perturb;
    }
}



    
 
// ************************************************************************* //
