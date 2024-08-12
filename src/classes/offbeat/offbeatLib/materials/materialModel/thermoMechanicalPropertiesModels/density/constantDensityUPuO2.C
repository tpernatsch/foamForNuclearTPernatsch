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

#include "constantDensityUPuO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantDensityUPuO2, 0);
    addToRunTimeSelectionTable
    (
        densityModel, 
        constantDensityUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantDensityUPuO2::constantDensityUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    densityModel(mesh, dict, defaultModel),       
    densityFrac_(dict.lookupOrDefault("densityFraction", 0.945)),
    theoreticalDensity_(dict.lookupOrDefault<scalar>("theoreticalDensity", 10430.0))
{
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantDensityUPuO2::~constantDensityUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::constantDensityUPuO2::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{
    // Compute initial porosity
    const scalar p0 = 1-densityFrac_;

    // Account for density change due to porosity migration if field exists
    const scalarField& p = mesh_.foundObject<volScalarField>("porosity")
    ? mesh_.lookupObject<volScalarField>("porosity").internalField().field()
    : scalarField(mesh_.C().size(), p0);

    forAll(addr, i)
    {   
        const label cellI = addr[i];

        const scalar currentDensityFrac(1-p[cellI]);

        sf[cellI] = currentDensityFrac * theoreticalDensity_;
    }


}

// ************************************************************************* //
