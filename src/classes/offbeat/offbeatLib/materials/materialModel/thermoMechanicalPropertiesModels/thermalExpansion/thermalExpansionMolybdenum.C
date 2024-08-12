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

#include "thermalExpansionMolybdenum.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionMolybdenum, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionMolybdenum, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionMolybdenum::thermalExpansionMolybdenum
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    par1(4.985E-06),
    par2(6.667E-10),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        par1 = thermalExpansionDict.lookupOrDefault<scalar>("par1", 4.985E-06);
        par2 = thermalExpansionDict.lookupOrDefault<scalar>("par2", 6.667E-10);  

        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);    
    }

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionMolybdenum::~thermalExpansionMolybdenum()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionMolybdenum::correct
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{      
    forAll(addr, i)
    {   
        const label cellI = addr[i];      

        const scalar Ti = T[cellI];

        symmTensor nominalValue((par1 + par2*Ti)*(Ti - 273.15)*I);

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
