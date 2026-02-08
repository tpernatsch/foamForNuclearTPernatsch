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

#include "heatCapacityBanerjee1515Ti.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityBanerjee1515Ti, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityBanerjee1515Ti, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityBanerjee1515Ti::heatCapacityBanerjee1515Ti
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, baseDict, defaultModel),
    par1(431),
    par2(0.177),
    par3(8.72e-5)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'heatCapacity' dictionary
    if((dict.found("heatCapacity")) || (dict.dictName() == "heatCapacity"))
    {
        const dictionary& heatCapacityDict = dict.found("heatCapacity")?
        dict.subDict("heatCapacity"):
        dict;

        par1 = heatCapacityDict.lookupOrDefault<scalar>("par1", 431);
        par2 = heatCapacityDict.lookupOrDefault<scalar>("par2", 0.177);
        par3 = heatCapacityDict.lookupOrDefault<scalar>("par3", 8.72e-5);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityBanerjee1515Ti::~heatCapacityBanerjee1515Ti()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacityBanerjee1515Ti::correct
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

        if (Ti < 293 || Ti > 1273)
        {
            WarningInFunction
                << "Foam::heatCapacityBanerjee1515Ti::correct() :\n"
                << "Supplied temperature, " << Ti 
                << ", is outside of range 293 < T < 1273 K " << endl;
        }

        const scalar nominalValue = par1 + par2 * Ti + par3 * pow(Ti,-2);

        sf[cellI] = nominalValue;

    }
}




    
 
// ************************************************************************* //
