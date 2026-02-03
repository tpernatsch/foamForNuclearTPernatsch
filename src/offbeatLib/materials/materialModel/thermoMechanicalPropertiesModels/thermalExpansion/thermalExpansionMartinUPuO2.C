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
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    par1_(9.828e-6),
    par2_(-6.39e-10),
    par3_(1.33e-12),
    par4_(-1.757e-17),
    par5_(1.1833e-5),
    par6_(-5.013e-9),
    par7_(3.756e-12),
    par8_(-6.125e-17),    
    OM_
    (
        dict.lookupOrDefault<scalar>
        (
            "oxygenMetalRatio",
            baseDict.lookupOrDefault<scalar>("oxygenMetalRatio", 2.0)
        )
    ),
    x_(2.0 - OM_)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'thermalExpansion' dictionary
    if((dict.found("thermalExpansion")) || (dict.dictName() == "thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.found("thermalExpansion")?
        dict.subDict("thermalExpansion"):
        dict;

        par1_ = thermalExpansionDict.lookupOrDefault<scalar>("par1", 9.828e-6);
        par2_ = thermalExpansionDict.lookupOrDefault<scalar>("par2", -6.39e-10);
        par3_ = thermalExpansionDict.lookupOrDefault<scalar>("par3", 1.33e-12);
        par4_ = thermalExpansionDict.lookupOrDefault<scalar>("par4", -1.757e-17);
        par5_ = thermalExpansionDict.lookupOrDefault<scalar>("par5", 1.1833e-5);
        par6_ = thermalExpansionDict.lookupOrDefault<scalar>("par6", -5.013e-9);
        par7_ = thermalExpansionDict.lookupOrDefault<scalar>("par7", 3.756e-12);
        par8_ = thermalExpansionDict.lookupOrDefault<scalar>("par8", -6.125e-17);
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
    // Get reference temperature
    const scalar Tref(Tref_);

    const scalar alphaRef = 
    (par1_ + par2_*Tref + par3_*pow(Tref,2) + par4_*pow(Tref,3)) * (1 + 3.98*x_);

    forAll(addr, i)
    {   
        const label cellI = addr[i];    

        const scalar Ti = T[cellI];

        const scalar alpha = ( Ti <= 923 )
        ? (par1_ + par2_*Ti + par3_*pow(Ti,2) + par4_*pow(Ti,3)) * (1 + 3.98*x_)
        : (par5_ + par6_*Ti + par7_*pow(Ti,2) + par8_*pow(Ti,3)) * (1 + 3.98*x_); 
        
        
        // Compute thermal expansion tensor
        // symmTensor nominalValue = alpha*(Ti-Tref)*I;
        symmTensor nominalValue = (alpha*Ti-alphaRef*Tref)*I;

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
