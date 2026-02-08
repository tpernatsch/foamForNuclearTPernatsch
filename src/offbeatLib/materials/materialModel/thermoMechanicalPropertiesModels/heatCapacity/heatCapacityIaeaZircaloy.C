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

#include "heatCapacityIaeaZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityIaeaZircaloy, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityIaeaZircaloy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityIaeaZircaloy::heatCapacityIaeaZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, baseDict, defaultModel),
    par1(255.66),
    par2(0.1024),
    par3(597.1),
    par4(0.4088),
    par5(1.565E-4),
    par6(1058.4),
    par7(1213.8),
    par8(2.0),
    par9(719.61)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'heatCapacity' dictionary
    if((dict.found("heatCapacity")) || (dict.dictName() == "heatCapacity"))
    {
        const dictionary& heatCapacityDict = dict.found("heatCapacity")?
        dict.subDict("heatCapacity"):
        dict;

        par1 = heatCapacityDict.lookupOrDefault<scalar>("par1", 255.66);
        par2 = heatCapacityDict.lookupOrDefault<scalar>("par2", 0.1024);
        par3 = heatCapacityDict.lookupOrDefault<scalar>("par3", 597.1);
        par4 = heatCapacityDict.lookupOrDefault<scalar>("par4", 0.4088);
        par5 = heatCapacityDict.lookupOrDefault<scalar>("par5", 1.565E-4);
        par6 = heatCapacityDict.lookupOrDefault<scalar>("par6", 1058.4);
        par7 = heatCapacityDict.lookupOrDefault<scalar>("par7", 1213.8);
        par8 = heatCapacityDict.lookupOrDefault<scalar>("par8", 2.0);
        par9 = heatCapacityDict.lookupOrDefault<scalar>("par9", 719.61);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityIaeaZircaloy::~heatCapacityIaeaZircaloy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::heatCapacityIaeaZircaloy::correct
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

        if (Ti < 273 || Ti > 2000)
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", is outside of range 273 < T < 2000 K";
        }

        scalar Cp1 = par1 + par2*Ti;
        scalar Cp2 = par3 - par4*Ti + par5*Ti*Ti;
        scalar f = par6*exp(-pow(Ti - par7,par8)/par9);

        scalar nominalValue(0.0);
        
        if (Ti < 1100)
        {
            nominalValue = Cp1;
        }
        else if (Ti < 1213.8)
        {
            nominalValue = Cp1 + f;
        }
        else if (Ti < 1320)
        {
            nominalValue = Cp2 + f;
        }
        else
        {
            nominalValue = Cp2;
        }

        sf[cellI] = nominalValue;

    }
}




    
 
// ************************************************************************* //
