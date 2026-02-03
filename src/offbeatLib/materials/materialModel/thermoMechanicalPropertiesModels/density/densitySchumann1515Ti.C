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

#include "densitySchumann1515Ti.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densitySchumann1515Ti, 0);
    addToRunTimeSelectionTable
    (
        densityModel, 
        densitySchumann1515Ti, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densitySchumann1515Ti::densitySchumann1515Ti
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    densityModel(mesh, dict, baseDict, defaultModel),       
    rho0_(7900.0)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'density' dictionary
    if((dict.found("density")) || (dict.dictName() == "density"))
    {
        const dictionary& densityDict = dict.found("density")?
        dict.subDict("density"):
        dict;

        rho0_ = densityDict.lookupOrDefault<scalar>("par1", 7900.0);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densitySchumann1515Ti::~densitySchumann1515Ti()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::densitySchumann1515Ti::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{
    scalar linExpCoeff;
    
    forAll(addr, i)
    {   
        const label cellI = addr[i];

        const scalar Ti = T[cellI] - 273.15;

        linExpCoeff = -3.101e-4 + 1.545e-5*Ti + 2.75e-9*pow(Ti,2);

        const scalar nominalValue =
        rho0_ * pow(1/(1+linExpCoeff),3);

        sf[cellI] = nominalValue;
    }
}

    

// ************************************************************************* //
