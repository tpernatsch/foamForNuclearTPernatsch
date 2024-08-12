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

#include "constantDensityUO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantDensityUO2, 0);
    addToRunTimeSelectionTable
    (
        densityModel, 
        constantDensityUO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantDensityUO2::constantDensityUO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    densityModel(mesh, dict, defaultModel),       
    densityFrac_(dict.lookupOrDefault("densityFraction", 0.95)),
    theoreticalDensity_(dict.lookupOrDefault<scalar>("theoreticalDensity", 10960.0))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantDensityUO2::~constantDensityUO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::constantDensityUO2::correct
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

        sf[cellI] = densityFrac_ * theoreticalDensity_;
    }
}

// ************************************************************************* //
