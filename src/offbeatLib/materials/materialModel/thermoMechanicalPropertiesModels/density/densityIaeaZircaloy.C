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

#include "densityIaeaZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densityIaeaZircaloy, 0);
    addToRunTimeSelectionTable
    (
        densityModel, 
        densityIaeaZircaloy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densityIaeaZircaloy::densityIaeaZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    densityModel(mesh, dict, baseDict, defaultModel),       
    par1(6595.2),
    par2(0.1477),
    par3(6690.0),
    par4(0.1855)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'density' dictionary
    if((dict.found("density")) || (dict.dictName() == "density"))
    {
        const dictionary& densityDict = dict.found("density")?
        dict.subDict("density"):
        dict;

        par1 = densityDict.lookupOrDefault<scalar>("par1", 6595.2);
        par2 = densityDict.lookupOrDefault<scalar>("par2", 0.1477);
        par3 = densityDict.lookupOrDefault<scalar>("par3", 6690.0);
        par4 = densityDict.lookupOrDefault<scalar>("par4", 0.1855);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densityIaeaZircaloy::~densityIaeaZircaloy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::densityIaeaZircaloy::correct
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
        
        scalar rho1 = par1 - par2*Ti;
        scalar rho2 = par3 - par4*Ti; 

        scalar nominalValue(0.0);
        
        if (Ti < 1083)
        {
            nominalValue = rho1;
        }
        else if (Ti < 1144)
        {
            nominalValue = rho1 + (rho2-rho1)*(Ti-1083)/(1144-1083);
        }
        else if (Ti < 1800)
        {
            nominalValue = rho2;
        }
        else
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", above maximum of 1800 K";
            
            nominalValue = 0;
        }

        sf[cellI] = nominalValue;
    }
}

    

// ************************************************************************* //
