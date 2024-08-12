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

#include "densityIAEAZy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densityIAEAZy, 0);
    addToRunTimeSelectionTable
    (
        densityModel, 
        densityIAEAZy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densityIAEAZy::densityIAEAZy
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    densityModel(mesh, dict, defaultModel),       
    par1(6595.2),
    par2(0.1477),
    par3(6690.0),
    par4(0.1855),
    perturb(1)
{
    if(dict.found("density"))
    {
        const dictionary& densityDict = dict.subDict("density");

        par1 = densityDict.lookupOrDefault<scalar>("par1", 6595.2);
        par2 = densityDict.lookupOrDefault<scalar>("par2", 0.1477);
        par3 = densityDict.lookupOrDefault<scalar>("par3", 6690.0);
        par4 = densityDict.lookupOrDefault<scalar>("par4", 0.1855);
        
        perturb = densityDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densityIAEAZy::~densityIAEAZy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::densityIAEAZy::correct
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

        sf[cellI] = nominalValue*perturb;
    }
}

    

// ************************************************************************* //
