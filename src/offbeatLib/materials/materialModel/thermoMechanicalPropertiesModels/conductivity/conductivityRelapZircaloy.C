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

#include "conductivityRelapZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityRelapZircaloy, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityRelapZircaloy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityRelapZircaloy::conductivityRelapZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, baseDict, defaultModel),     
    par1(7.51),
    par2(2.09e-2),
    par3(1.45e-5),
    par4(7.67e-9)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'conductivity' dictionary
    if((dict.found("conductivity")) || (dict.dictName() == "conductivity"))
    {
        const dictionary& conductivityDict = dict.found("conductivity")?
        dict.subDict("conductivity"):
        dict;

        par1 = conductivityDict.lookupOrDefault<scalar>("par1", 7.51);
        par2 = conductivityDict.lookupOrDefault<scalar>("par2", 2.09e-2);
        par3 = conductivityDict.lookupOrDefault<scalar>("par3", 1.45e-5);
        par4 = conductivityDict.lookupOrDefault<scalar>("par4", 7.67e-9);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityRelapZircaloy::~conductivityRelapZircaloy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityRelapZircaloy::correct
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

        //- The validity range is changed by one degree on each side to allow 
        //- for rounding errors
        if (Ti < 290 || Ti > 1801)
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", out of range 290 < T < 1800 K";
        }
       
        const scalar nominalValue =  par1 + par2*Ti - par3*Ti*Ti + par4*Ti*Ti*Ti;

        sf[cellI] = nominalValue;
    }
}


    


// ************************************************************************* //
