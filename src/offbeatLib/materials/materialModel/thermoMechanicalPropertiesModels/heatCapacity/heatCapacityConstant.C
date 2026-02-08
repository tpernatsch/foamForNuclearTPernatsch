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

#include "heatCapacityConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatCapacityConstant, 0);
    addToRunTimeSelectionTable
    (
        heatCapacityModel, 
        heatCapacityConstant, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatCapacityConstant::heatCapacityConstant
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    heatCapacityModel(mesh, dict, baseDict, defaultModel),
    cp_()
{
    if(dict.dictName() == "heatCapacity")
    {
        dict.lookup("value") >> cp_;
    }
    else
    {
        // For retrocompatibility
        cp_ = dimensionedScalar(dict.lookup("Cp")).value();
        
        WarningIn("heatCapacityConstant::New(const fvMesh&, const dictionary&, const word&)")
        << "Keyword 'Cp' as dimensionedScalar in the material dictionary is deprecated for constant heatCapacity model. " << nl
        << "Please create a 'heatCapacity' sub-dictionary with 'constant' as type, and use keyword 'value' to assign the constant heatCapacity value." << endl;
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::heatCapacityConstant::~heatCapacityConstant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatCapacityConstant::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{
    forAll(addr, i)
    {
        sf[addr[i]] = cp_;
    }
}




    
 
// ************************************************************************* //
