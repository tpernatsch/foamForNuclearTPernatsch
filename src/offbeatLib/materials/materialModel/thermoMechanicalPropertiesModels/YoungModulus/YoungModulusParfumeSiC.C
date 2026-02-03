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

#include "YoungModulusParfumeSiC.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusParfumeSiC, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel,
        YoungModulusParfumeSiC,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusParfumeSiC::YoungModulusParfumeSiC
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, baseDict, defaultModel),
    temp1(25.0),
    temp2(940.0),
    temp3(1215.0),
    temp4(1600.0),
    E1(428.0),
    E2(375.0),
    E3(340.0),
    E4(198.0)
{

    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'YoungModulus' dictionary
    if((dict.found("YoungModulus")) || (dict.dictName() == "YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.found("YoungModulus")?
        dict.subDict("YoungModulus"):
        dict;

        temp1 = YoungModulusDict.lookupOrDefault<scalar>("temp1", 25.0);
        temp2 = YoungModulusDict.lookupOrDefault<scalar>("temp2", 940.0);
        temp3 = YoungModulusDict.lookupOrDefault<scalar>("temp3", 1215.0);
        temp4 = YoungModulusDict.lookupOrDefault<scalar>("temp4", 1600.0);
        E1 = YoungModulusDict.lookupOrDefault<scalar>("E1", 428.0);
        E2 = YoungModulusDict.lookupOrDefault<scalar>("E2", 375.0);
        E3 = YoungModulusDict.lookupOrDefault<scalar>("E3", 340.0);
        E4 = YoungModulusDict.lookupOrDefault<scalar>("E4", 198.0);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusParfumeSiC::~YoungModulusParfumeSiC()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusParfumeSiC::correct
(
    scalarField& sf,
    const scalarField& T,
    const labelList& addr
)
{

    forAll(addr, i)
    {
        const label cellI = addr[i];

        const scalar Ti = T[cellI] - 273.15;  //-temperature in Celsius

        scalar nominalValue(0.0);

        //- Interpolation
        if (Ti <= temp1)
        {
          nominalValue = E1;
        }
        else if (Ti > temp1 && Ti <= temp2)
        {
          nominalValue = (E2-E1)/(temp2-temp1)*(Ti-temp1) + E1;
        }
        else if (Ti > temp2 && Ti <= temp3)
        {
          nominalValue = (E3-E2)/(temp3-temp2)*(Ti-temp2) + E2;
        }
        else if (Ti > temp3 && Ti <= temp4)
        {
          nominalValue = (E4-E3)/(temp4-temp3)*(Ti-temp3) + E3;
        }
        else
        {
          nominalValue = E4;
        }

        sf[cellI] = 1e9*nominalValue;
    }
}




// ************************************************************************* //
