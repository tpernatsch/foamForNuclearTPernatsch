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

#include "thermalExpansionMartinUPuO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionMartinUPuO2, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionMartinUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionMartinUPuO2::thermalExpansionMartinUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict, 
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, defaultModel),
    par1_(-2.66e-3),
    par2_(9.802e-6),
    par3_(-2.705e-10),
    par4_(4.391e-13),
    par5_(-3.28e-3),
    par6_(1.179e-5),
    par7_(-2.429e-9),
    par8_(1.219e-12),
    perturb(1.0)
{
    if(dict.found("thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.subDict("thermalExpansion");

        par1_ = thermalExpansionDict.lookupOrDefault<scalar>("par1_", -2.66e-3);
        par2_ = thermalExpansionDict.lookupOrDefault<scalar>("par2_", 9.802e-6);
        par3_ = thermalExpansionDict.lookupOrDefault<scalar>("par3_", -2.705e-10);
        par4_ = thermalExpansionDict.lookupOrDefault<scalar>("par4_", 4.391e-13);
        par5_ = thermalExpansionDict.lookupOrDefault<scalar>("par5_", -3.28e-3);
        par6_ = thermalExpansionDict.lookupOrDefault<scalar>("par6_", 1.179e-5);
        par7_ = thermalExpansionDict.lookupOrDefault<scalar>("par7_", -2.429e-9);
        par8_ = thermalExpansionDict.lookupOrDefault<scalar>("par8_", 1.219e-12);

        perturb = thermalExpansionDict.lookupOrDefault<scalar>("perturb", 1.0);   
        
    }
}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionMartinUPuO2::~thermalExpansionMartinUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionMartinUPuO2::correct
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

        symmTensor nominalValue = ( Ti<923 )
            ? (par1_ + par2_*Ti + par3_*pow(Ti,2) + par4_*pow(Ti,3))*I
            : (par5_ + par6_*Ti + par7_*pow(Ti,2) + par8_*pow(Ti,3))*I; 
        

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
