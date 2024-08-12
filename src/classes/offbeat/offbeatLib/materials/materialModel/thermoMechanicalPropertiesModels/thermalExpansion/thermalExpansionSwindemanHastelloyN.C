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

#include "thermalExpansionSwindemanHastelloyN.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionSwindemanHastelloyN, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel,
        thermalExpansionSwindemanHastelloyN,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionSwindemanHastelloyN::thermalExpansionSwindemanHastelloyN
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    par1(0.005291),
    par2(9.682),
    par3(107.8),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        par1 = thermalExpansionDict.lookupOrDefault<scalar>("par1", 0.005291);
        par2 = thermalExpansionDict.lookupOrDefault<scalar>("par2", 9.682);
        par3 = thermalExpansionDict.lookupOrDefault<scalar>("par3", 107.8);
 
        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);
    }

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionSwindemanHastelloyN::~thermalExpansionSwindemanHastelloyN()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionSwindemanHastelloyN::correct
(
 symmTensorField& sf,
 const scalarField& T,
 const labelList& addr
)
const
{
    const scalar Tref = Tref_.value();
    const symmTensor referenceValue =
    (par1*(Tref-273.15)*(Tref-273.15) + par2*(Tref-273.15) + par3)*(1.0e-6)*I;

    forAll(addr, i)
    {
        const label cellI = addr[i];
        const scalar Ti = T[cellI];

        symmTensor nominalValue =
        (par1*(Ti-273.15)*(Ti-273.15) + par2*(Ti-273.15) + par3)*(1.0e-6)*I
        - referenceValue;
     
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

        sf[cellI] = nominalValue*perturb;
    }
}


// ************************************************************************* //