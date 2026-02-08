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
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    par1_(-3.9735e-4),
    par2_( 8.4955e-6),
    par3_( 2.15130e-9),
    par4_( 3.7143e-16),
    par5_(-4.972e-4),
    par6_( 7.107e-6),
    par7_( 2.581e-9),
    par8_( 1.140e-13),
    N_U_(mesh_.lookupObject<volScalarField>("N_U")),
    N_Pu_(mesh_.lookupObject<volScalarField>("N_Pu")),
    N_Am_(nullptr),
    N_Np_(nullptr)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'thermalExpansion' dictionary
    if((dict.found("thermalExpansion")) || (dict.dictName() == "thermalExpansion"))
    {
        const dictionary& thermalExpansionDict = dict.found("thermalExpansion")?
        dict.subDict("thermalExpansion"):
        dict;

        par1_ = thermalExpansionDict.lookupOrDefault<scalar>("par1", -3.9735e-4);
        par2_ = thermalExpansionDict.lookupOrDefault<scalar>("par2", 8.4955e-6);
        par3_ = thermalExpansionDict.lookupOrDefault<scalar>("par3", 2.15130e-9);
        par4_ = thermalExpansionDict.lookupOrDefault<scalar>("par4", 3.7143e-16);
        par5_ = thermalExpansionDict.lookupOrDefault<scalar>("par5", -4.972e-4);
        par6_ = thermalExpansionDict.lookupOrDefault<scalar>("par6",  7.107e-6);
        par7_ = thermalExpansionDict.lookupOrDefault<scalar>("par7",  2.581e-9);
        par8_ = thermalExpansionDict.lookupOrDefault<scalar>("par8",  1.140e-13);
    }

    // Check if minor actinides have been initialized by materials
    if( mesh_.foundObject<volScalarField>("N_Am") )
    {
        N_Am_ = &mesh_.lookupObjectRef<volScalarField>("N_Am");
    }
    if( mesh_.foundObject<volScalarField>("N_Np") )
    {
        N_Np_ = &mesh_.lookupObjectRef<volScalarField>("N_Np");
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

        // Temperature in C
        const scalar Ti = T[cellI] - 273.15;

        // Number densities for heavy metals in fuel
        const scalar N_U = N_U_[cellI];
        const scalar N_Pu = N_Pu_[cellI];
        const scalar N_Am = (N_Am_ == nullptr) ? 0.0 : N_Am_->internalField()[cellI];
        const scalar N_Np = (N_Np_ == nullptr) ? 0.0 : N_Np_->internalField()[cellI];

        // Compute concentrations atoms/atomsHM
        const scalar conc_at_Pu = N_Pu / (N_Pu + N_Am + N_Np + N_U);    

        // Conversion factors from atoms/atomsHM -> mass/mass_fuel
        // Approximate conversion factors are used, computed in this way:
        // 
        // For Pu:
        // at.Pu/at.HM = k_Pu * m_Pu/m_HM
        // Assuming: 
        // MM_Pu ~ 239 g/mol, OM_ratio = 2, MM_fuel ~ 238.5 g/mol, MM_O = 16:
        // => k_Pu = (MM_HM + OM_ratio * MM_O) / MM_Pu
        const scalar k_Pu = 1.13;

        // Assign concentrations in mass/mass_fuel
        const scalar c_w_Pu = conc_at_Pu / k_Pu;

        // Linear thermal expansion of PuO2
        scalar dL_L_PuO2 = 
            (par1_ + par2_*Ti + par3_*pow(Ti,2) + par4_*pow(Ti,3));

        // Linear thermal expansion of UO2
        scalar dL_L_UO2 = 
            (par5_ + par6_*Ti + par7_*pow(Ti,2) + par8_*pow(Ti,3));
        
        // Weigthed average
        scalar dL_L = c_w_Pu * dL_L_PuO2 + ( 1-c_w_Pu) * dL_L_UO2;

        // Assign thermal expansion
        symmTensor nominalValue = dL_L * I;

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
