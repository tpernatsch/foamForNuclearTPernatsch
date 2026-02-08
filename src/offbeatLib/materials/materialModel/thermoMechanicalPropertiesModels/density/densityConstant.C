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

#include "densityConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(densityConstant, 0);
    addToRunTimeSelectionTable
    (
        densityModel, 
        densityConstant, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::densityConstant::densityConstant
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    densityModel(mesh, dict, baseDict, defaultModel),
    rho_()
{
    if(dict.dictName() == "density")
    {
        dict.lookup("value") >> rho_;
    }
    else
    {
        // For retrocompatibility
        rho_ = dimensionedScalar(dict.lookup("rho")).value();
        
        WarningIn("densityConstant::New(const fvMesh&, const dictionary&, const word&)")
        << "Keyword 'rho' as dimensionedScalar in the material dictionary is deprecated for constant density model. " << nl
        << "Please create a 'density' sub-dictionary with 'constant' as type, and use keyword 'value' to assign the constant density value." << endl;    
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::densityConstant::~densityConstant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::densityConstant::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{
    forAll(addr, i)
    {
        sf[addr[i]] = rho_;
    }
}

    

// ************************************************************************* //
