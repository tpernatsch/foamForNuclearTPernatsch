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

#include "thermalExpansionMatproUPuO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionMatproUPuO2, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionMatproUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionMatproUPuO2::thermalExpansionMatproUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    par1_(-3.9735e-4),
    par2_(8.4955e-6),
    par3_(2.15130e-9),
    par4_(3.7143e-16),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        par1_ = thermalExpansionDict.lookupOrDefault<scalar>("par1", -3.9735e-4);
        par1_ = thermalExpansionDict.lookupOrDefault<scalar>("par2", 8.4955e-6);
        par1_ = thermalExpansionDict.lookupOrDefault<scalar>("par3", 2.15130e-9);
        par1_ = thermalExpansionDict.lookupOrDefault<scalar>("par4", 3.7143e-16);

        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);   
        
    }
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionMatproUPuO2::~thermalExpansionMatproUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionMatproUPuO2::correct
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
        
        const symmTensor nominalValue = 
        (par1_ + par2_*Ti + par3_*pow(Ti,2) + par4_*pow(Ti,3))*I;

        sf[cellI] = nominalValue*perturb;

    }
}   
 
// ************************************************************************* //
