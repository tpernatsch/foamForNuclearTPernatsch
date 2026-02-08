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

#include "YoungModulusWatrousHastelloyN.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusWatrousHastelloyN, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel,
        YoungModulusWatrousHastelloyN,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusWatrousHastelloyN::YoungModulusWatrousHastelloyN
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, baseDict, defaultModel),
    par1(9.944e-08),
    par2(1.178e-4),
    par3(1.033e-1),
    par4(220.9)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'YoungModulus' dictionary
    if((dict.found("YoungModulus")) || (dict.dictName() == "YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.found("YoungModulus")?
        dict.subDict("YoungModulus"):
        dict;

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 9.944e-08);
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 1.178e-4);
        par3 = YoungModulusDict.lookupOrDefault<scalar>("par3", 1.033e-1);
        par4 = YoungModulusDict.lookupOrDefault<scalar>("par4", 220.9);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusWatrousHastelloyN::~YoungModulusWatrousHastelloyN()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusWatrousHastelloyN::correct
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


        scalar nominalValue(0.0);
        
        // Below 1000C
        if (Ti < 1273.15)
        {
            scalar K1 = -par1*(Ti-273.15)*(Ti-273.15)*(Ti-273.15);
            scalar K2 = par2*(Ti-273.15)*(Ti-273.15);
            scalar K3 = -par3*(Ti-273.15) + par4;
            
            nominalValue = (K1 + K2 + K3)*1e9;
        }
    
        else
        
        {
                WarningInFunction
                    << "Supplied temperature, " << Ti << ", is outside of range: T < 1273.15 K";
        }

        sf[cellI] = nominalValue;
    }
}