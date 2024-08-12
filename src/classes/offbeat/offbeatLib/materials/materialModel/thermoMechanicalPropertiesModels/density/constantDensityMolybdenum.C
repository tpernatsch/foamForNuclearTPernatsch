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

#include "constantDensityMolybdenum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantDensityMolybdenum, 0);
    addToRunTimeSelectionTable
    (
        densityModel, 
        constantDensityMolybdenum, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantDensityMolybdenum::constantDensityMolybdenum
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    densityModel(mesh, dict, defaultModel),       
    densityValue_(10280)
{
    if(dict.found("density"))
    {
        const dictionary& densityDict = dict.subDict("density");

        densityValue_ = 
        densityDict.lookupOrDefault<scalar>("densityValue", 10280);  
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantDensityMolybdenum::~constantDensityMolybdenum()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::constantDensityMolybdenum::correct
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
        
        sf[cellI] =  densityValue_;
    }
}
// ************************************************************************* //
