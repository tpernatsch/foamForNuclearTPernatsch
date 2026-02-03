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

#include "YoungModulusHofmanD9.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusHofmanD9, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel, 
        YoungModulusHofmanD9, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusHofmanD9::YoungModulusHofmanD9
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, baseDict, defaultModel),
    par1(2.01e5),
    par2(79.29)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'YoungModulus' dictionary
    if((dict.found("YoungModulus")) || (dict.dictName() == "YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.found("YoungModulus")?
        dict.subDict("YoungModulus"):
        dict;

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 2.01e5);        
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 79.29);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusHofmanD9::~YoungModulusHofmanD9()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusHofmanD9::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{ 

    forAll(addr, i)
    {   
        const label cellI = addr[i];        

        if (T[cellI] < 293 || T[cellI] > 1273)
        {
            WarningInFunction
                << "Foam::YoungModulusHofmanD9::correct() :\n"
                << "Supplied temperature, " << T[cellI] 
                << ", is outside of range 293 < T < 1273 K " << endl;
        }

        const scalar Ti = T[cellI] - 273.15;  

        const scalar nominalValue = (par1-par2*Ti)*1e6;
        
        sf[cellI] = nominalValue;
    }
}




// ************************************************************************* //
