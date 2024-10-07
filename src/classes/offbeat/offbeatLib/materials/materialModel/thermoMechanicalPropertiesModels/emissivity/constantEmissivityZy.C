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

#include "constantEmissivityZy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constantEmissivityZy, 0);
    addToRunTimeSelectionTable
    (
        emissivityModel, 
        constantEmissivityZy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantEmissivityZy::constantEmissivityZy
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    emissivityModel(mesh, dict, defaultModel),
    emissivityValue_(0.808642)
{
    if(dict.found("emissivity"))
    {
        const dictionary& emissivityDict = dict.subDict("emissivity");

        emissivityValue_ = 
        emissivityDict.lookupOrDefault<scalar>("emissivityValue", 0.808642);        
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constantEmissivityZy::~constantEmissivityZy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::constantEmissivityZy::correct
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

        if (Ti > 1500)
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", out of range T < 1500 K";
        }
        
        sf[cellI] = emissivityValue_;
    }
}




    
 
// ************************************************************************* //
