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

#include "thermalExpansionLemehovUPuO2.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermalExpansionLemehovUPuO2, 0);
    addToRunTimeSelectionTable
    (
        thermalExpansionModel, 
        thermalExpansionLemehovUPuO2, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermalExpansionLemehovUPuO2::thermalExpansionLemehovUPuO2
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    thermalExpansionModel(mesh, dict, baseDict, defaultModel),
    b0_(-0.3080),
    b1_(3.4303),
    b2_(-1.9157),
    b3_(3.4636),
    by_(3.98),
    OM_(mesh_.lookupObject<volScalarField>("oxygenMetalRatio")),
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

        b0_ = thermalExpansionDict.lookupOrDefault<scalar>("b0", -0.3080);
        b1_ = thermalExpansionDict.lookupOrDefault<scalar>("b1",  3.4303);
        b2_ = thermalExpansionDict.lookupOrDefault<scalar>("b2", -1.9157);
        b3_ = thermalExpansionDict.lookupOrDefault<scalar>("b3",  3.4636);
        by_ = thermalExpansionDict.lookupOrDefault<scalar>("by",  3.98);
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

Foam::thermalExpansionLemehovUPuO2::~thermalExpansionLemehovUPuO2()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::thermalExpansionLemehovUPuO2::correct
(
    symmTensorField& sf, 
    const scalarField& T, 
    const labelList& addr
)
const
{
    const volScalarField& Bu(mesh_.lookupObject<volScalarField>("Bu"));

    forAll(addr, i)
    {
        const label cellI = addr[i];    

        const scalar Ti = T[cellI];

        const scalar BuGWd_tHM = Bu.internalField()[cellI] / 1000;

        // Number densities for heavy metals in fuel
        const scalar N_U = N_U_[cellI];
        const scalar N_Pu = N_Pu_[cellI];
        const scalar N_Am = (N_Am_ == nullptr) ? 0.0 : N_Am_->internalField()[cellI];
        const scalar N_Np = (N_Np_ == nullptr) ? 0.0 : N_Np_->internalField()[cellI];
        
        // Compute concentrations atoms/atomsHM
        const scalar conc_at_Pu = N_Pu / (N_Pu + N_Am + N_Np + N_U);    
        const scalar conc_at_Am = N_Am / (N_Pu + N_Am + N_Np + N_U);

        // Conversion factors from atoms/atomsHM -> mass/mass_fuel
        // Approximate conversion factors are used, computed in this way:
        // 
        // For Pu:
        // at.Pu/at.HM = k_Pu * m_Pu/m_HM
        // Assuming: 
        // MM_Pu ~ 239 g/mol, OM_ratio = 2, MM_fuel ~ 238.5 g/mol, MM_O = 16:
        // => k_Pu = (MM_HM + OM_ratio * MM_O) / MM_Pu
        const scalar k_Pu = 1.13;
        const scalar k_Am = 1.12;

        // Assign concentrations in mass/mass_fuel
        const scalar c_w_Pu = conc_at_Pu / k_Pu;
        const scalar c_w_Am = conc_at_Am / k_Am;

        // // Compute melting temperature according to KONNO (https://doi.org/10.1080/18811248.1998.9733897)
        // const scalar Tm = 3133.8 
        //             - 460  * c_w_Pu 
        //             - 980  * c_w_Am 
        //             - 0.66 * BuGWd_tHM + 0.0013*pow(BuGWd_tHM,2.0);

        // Compute melting temperature according to MAGNI (https://doi.org/10.1016/j.jnucmat.2021.153312)
        const scalar Tm0 = 3147
                    -  364.85  * c_w_Pu 
                    -  1014.15 * (2.0-OM_[cellI])
                    -  329.5   * c_w_Am;
        
        const scalar Tm = 2964.94 + (Tm0 - 2964.94) * exp(-BuGWd_tHM/24.25);
        
        // Compute dL/L0 (convert from % to -)
        const scalar dL_L0 = 0.01 * (b0_
                              + b1_*Ti/Tm
                              + b2_*pow(Ti/Tm,2.0)
                              + b3_*pow(Ti/Tm,3.0));

        // Compute dL/Ly accounting for deviation from soichiometry
        const scalar dL_Ly = dL_L0 * (1+by_*(2-OM_[cellI]));
        
        // Compute thermal expansion tensor
        symmTensor nominalValue = (dL_Ly) * I;

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
