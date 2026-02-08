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

#include "thermalExpansionConstant.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionConstant, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionConstant, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionConstant::thermalExpansionConstant
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    alpha_()
{
    if(dict.dictName() == "thermalExpansion")
    {
        dict.lookup("value") >> alpha_;
    }
    else
    {
        // For retrocompatibility
        alpha_ = dimensionedScalar(dict.lookup("alpha")).value();

        WarningIn("thermalExpansionConstant::New(const fvMesh&, const dictionary&, const word&)")
        << "Keyword 'alpha' as dimensionedScalar in the material dictionary is deprecated for constant thermalExpansion model. " << nl
        << "Please create a 'thermalExpansion' sub-dictionary with 'constant' as type, and use keyword 'value' to assign the constant thermalExpansion value." << endl;
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionConstant::~thermalExpansionConstant()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionConstant::correct
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{   
    forAll(addr, i)
    {
        symmTensor nominalValue(I*(alpha_*(T[addr[i]]-Tref_)));

        // NOTE: avoids instabilities when trying to simulate a material at 
        // constant temperature equal to Tref
        if
        (
            mag(nominalValue.xx()) < 1e-7 &&
            mag(nominalValue.yy()) < 1e-7 &&
            mag(nominalValue.zz()) < 1e-7
        )
        {
            nominalValue *= 0;
        }
        
        sf[addr[i]] = nominalValue;
    }
}

// ************************************************************************* //
