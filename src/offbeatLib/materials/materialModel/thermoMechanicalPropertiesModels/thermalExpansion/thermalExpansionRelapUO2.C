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

#include "thermalExpansionRelapUO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionRelapUO2, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionRelapUO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionRelapUO2::thermalExpansionRelapUO2
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    K1(9.8e-6),
    K2(2.61e-3),
    K3(3.16e-1),
    ED(1.32e-19),
    k(1.38e-23)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'thermalExpansion' dictionary
    if((dict.found("thermalExpansion")) || (dict.dictName() == "thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.found("thermalExpansion")?
        dict.subDict("thermalExpansion"):
        dict;

        K1 = thermalExpansionDict.lookupOrDefault<scalar>("K1", 9.8e-6);
        K2 = thermalExpansionDict.lookupOrDefault<scalar>("K2", 2.61e-3);
        K3 = thermalExpansionDict.lookupOrDefault<scalar>("K3", 3.16e-1);
        ED = thermalExpansionDict.lookupOrDefault<scalar>("ED", 1.32e-19);
        k = thermalExpansionDict.lookupOrDefault<scalar>("k", 1.38e-23);        
    }
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionRelapUO2::~thermalExpansionRelapUO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionRelapUO2::correct
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{  
    const scalar Tref = Tref_;
    const symmTensor referenceValue = (K1*Tref - K2 + K3*exp(-ED/k/Tref))*I;

    forAll(addr, i)
    {   
        const label cellI = addr[i];    

        const scalar Ti = T[cellI];
        
        symmTensor nominalValue = 
        (K1*Ti - K2 + K3*exp(-ED/k/Ti))*I - referenceValue;

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

        sf[cellI] = nominalValue; 

    }
}   
 
// ************************************************************************* //
