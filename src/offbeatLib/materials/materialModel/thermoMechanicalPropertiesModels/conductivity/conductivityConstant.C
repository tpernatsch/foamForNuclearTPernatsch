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

#include "conductivityConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityConstant, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityConstant, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityConstant::conductivityConstant
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, baseDict, defaultModel),
    k_()
{
    if(dict.dictName() == "conductivity")
    {
        dict.lookup("value") >> k_;
    }
    else
    {
        // For retrocompatibility
        k_ = dimensionedScalar(dict.lookup("k")).value();
        
        WarningIn("conductivityConstant::New(const fvMesh&, const dictionary&, const word&)")
            << "Keyword 'k' as dimensionedScalar in the material dictionary is deprecated for constant conductivity model. " << nl
            << "Please create a 'conductivity' sub-dictionary with 'constant' as type, and use keyword 'value' to assign the constant conductivity value." << endl;
            
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::conductivityConstant::~conductivityConstant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::conductivityConstant::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{   
    forAll(addr, i)
    {
        sf[addr[i]] = k_;
    }
}

// ************************************************************************* //
