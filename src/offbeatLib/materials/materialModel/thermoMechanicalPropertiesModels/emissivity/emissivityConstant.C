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

#include "emissivityConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(emissivityConstant, 0);
    addToRunTimeSelectionTable
    (
        emissivityModel, 
        emissivityConstant, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::emissivityConstant::emissivityConstant
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    emissivityModel(mesh, dict, baseDict, defaultModel),
    emissivity_()
{
    if(dict.dictName() == "emissivity")
    {
        dict.lookup("value") >> emissivity_;
    }
    else
    {
        // For retrocompatibility
        emissivity_ = dimensionedScalar(dict.lookup("emissivity")).value();
        
        WarningIn("emissivityConstant::New(const fvMesh&, const dictionary&, const word&)")
        << "Keyword 'emissivity' as dimensionedScalar in the material dictionary is deprecated for constant emissivity model. " << nl
        << "Please create a 'emissivity' sub-dictionary with 'constant' as type, and use keyword 'value' to assign the constant emissivity value." << endl;
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::emissivityConstant::~emissivityConstant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::emissivityConstant::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{  
    forAll(addr, i)
    {
        sf[addr[i]] = emissivity_;
    }
}

    
 
// ************************************************************************* //
