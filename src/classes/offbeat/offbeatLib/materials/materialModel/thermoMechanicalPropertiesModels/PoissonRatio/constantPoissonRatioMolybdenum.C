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

#include "constantPoissonRatioMolybdenum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantPoissonRatioMolybdenum, 0);
    addToRunTimeSelectionTable
    (
        PoissonRatioModel, 
        constantPoissonRatioMolybdenum, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantPoissonRatioMolybdenum::constantPoissonRatioMolybdenum
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    PoissonRatioModel(mesh, dict, defaultModel),       
    PoissonRatioValue_(0.31)
{
    if(dict.found("PoissonRatio"))
    {
        const dictionary& PoissonRatioDict = dict.subDict("PoissonRatio");

        PoissonRatioValue_ = 
        PoissonRatioDict.lookupOrDefault<scalar>("PoissonRatioValue", 0.316);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantPoissonRatioMolybdenum::~constantPoissonRatioMolybdenum()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::constantPoissonRatioMolybdenum::correct
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
