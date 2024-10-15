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

#include "constantPoissonRatioZy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantPoissonRatioZy, 0);
    addToRunTimeSelectionTable
    (
        PoissonRatioModel, 
        constantPoissonRatioZy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantPoissonRatioZy::constantPoissonRatioZy
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    PoissonRatioModel(mesh, dict, defaultModel),       
    PoissonRatioValue_(0.3)
{
    if(dict.found("PoissonRatio"))
    {
        const dictionary& PoissonRatioDict = dict.subDict("PoissonRatio");

        PoissonRatioValue_ = 
        PoissonRatioDict.lookupOrDefault<scalar>("PoissonRatioValue", 0.3);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantPoissonRatioZy::~constantPoissonRatioZy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::constantPoissonRatioZy::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{    
    forAll(addr, i)
    {   
        const label cellI = addr[i];    
        
        sf[cellI] =  PoissonRatioValue_;
    }
}
// ************************************************************************* //
