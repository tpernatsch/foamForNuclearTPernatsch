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

#include "thermalExpansionMatproZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionMatproZircaloy, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionMatproZircaloy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionMatproZircaloy::thermalExpansionMatproZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    par1(4.441e-6),
    par2(1.238e-3),
    par3(6.721e-6),
    par4(2.073e-3),
    par5(9.7e-6),
    par6(1.10e-2),
    par7(9.4e-3)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'thermalExpansion' dictionary
    if((dict.found("thermalExpansion")) || (dict.dictName() == "thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.found("thermalExpansion")?
        dict.subDict("thermalExpansion"):
        dict;

        par1 = thermalExpansionDict.lookupOrDefault<scalar>("par1", 4.441e-6);
        par2 = thermalExpansionDict.lookupOrDefault<scalar>("par2", 1.238e-3);
        par3 = thermalExpansionDict.lookupOrDefault<scalar>("par3", 6.721e-6);
        par4 = thermalExpansionDict.lookupOrDefault<scalar>("par4", 2.073e-3);
        par5 = thermalExpansionDict.lookupOrDefault<scalar>("par5", 9.7e-6);
        par6 = thermalExpansionDict.lookupOrDefault<scalar>("par6", 1.10e-2);
        par7 = thermalExpansionDict.lookupOrDefault<scalar>("par7", 9.45e-3);    
    }

}
// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermalExpansionMatproZircaloy::~thermalExpansionMatproZircaloy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionMatproZircaloy::correct
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{    
    const scalar Tref = Tref_;
    const symmTensor alphaTref = setAlphaT(Tref);

    forAll(addr, i)
    {   
        const label cellI = addr[i];      

        const scalar Ti = T[cellI];

        symmTensor nominalValue(setAlphaT(Ti) - alphaTref);

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
 

const Foam::symmTensor Foam::thermalExpansionMatproZircaloy::setAlphaT(const scalar& Ti) const 
{
    Foam::symmTensor alphaT = I;

    // The validity range is changed by one degree on each side to allow 
    // for rounding errors
    // NOTE: lower bound from literature (290) extended to 273 for allowing 
    // thermal contraction below 290 K
    if (Ti < 272.9 || Ti > 1801)
    {        
        WarningInFunction
            << "Supplied temperature, " << Ti << ", is out of range 273 < T < 1800 K";
    }

    if (Ti < 1073) // Alpha phase
    {
        // // Axial expansion 
        // // NOTE we assume isotrpic behavior 
        // alphaT.zz() = par1*Ti - par2;  
        // Axial expansion
        alphaT.zz() = par3*Ti - par4;  
        alphaT.xx() = par3*Ti - par4;  // Radial expansion
        alphaT.yy() = par3*Ti - par4;  // Radial expansion
    }
    else if (Ti < 1273) // Interpolate between alpha and beta values
    {
        scalar TAlpha = 1073;
        scalar TBeta  = 1273;
        

        // // NOTE we assume isotrpic behavior 
        // scalar axAlpha = par1*TAlpha - par2;
        // scalar axBeta  = par5*TBeta - par6;
        scalar axAlpha = par3*TAlpha - par4;
        scalar axBeta  = par5*TBeta - par7;

        scalar radAlpha = par3*TAlpha - par4;
        scalar radBeta  = par5*TBeta - par7;
        
        
        alphaT.zz() = axAlpha*(TBeta - Ti)/(TBeta - TAlpha) + axBeta*(Ti - TAlpha)/(TBeta - TAlpha);
        alphaT.xx() = radAlpha*(TBeta - Ti)/(TBeta - TAlpha) + radBeta*(Ti - TAlpha)/(TBeta - TAlpha);
        alphaT.yy() = radAlpha*(TBeta - Ti)/(TBeta - TAlpha) + radBeta*(Ti - TAlpha)/(TBeta - TAlpha);
    }
    else // Beta phase
    {  
        // // Axial expansion 
        // // NOTE we assume isotrpic behavior 
        // alphaT.zz() = par5*Ti - par6;
        alphaT.zz() = par5*Ti - par7;
        alphaT.xx() = par5*Ti - par7;  // Radial expansion
        alphaT.yy() = par5*Ti - par7;  // Radial expansion
    } 

    return alphaT;
}

// ************************************************************************* //
