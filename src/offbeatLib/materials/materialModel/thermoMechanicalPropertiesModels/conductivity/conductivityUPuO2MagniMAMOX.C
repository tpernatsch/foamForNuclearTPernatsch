
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

#include "conductivityUPuO2MagniMAMOX.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUPuO2MagniMAMOX, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUPuO2MagniMAMOX, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUPuO2MagniMAMOX::conductivityUPuO2MagniMAMOX
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, baseDict, defaultModel),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    OM_(mesh_.lookupObject<volScalarField>("oxygenMetalRatio")),
    N_U_(mesh_.lookupObject<volScalarField>("N_U")),
    N_Pu_(mesh_.lookupObject<volScalarField>("N_Pu")),
    N_Am_(nullptr),
    N_Np_(nullptr),
    porosity_(mesh.lookupObject<volScalarField>("porosity")),
    kappaInf_(1.755),     
    phi_(128.75e3),     
    A0_(0.01926),     
    Ax_(1.06e-6),     
    APu_(2.63e-8),     
    AAm_(0.596),     
    ANp_(2.22e-14),     
    B0_(2.39e-4),     
    BPu_(1.37e-13),     
    BAm_(5.47e-4),     
    BNp_(2.48e-14),     
    D_(5.27e9),     
    E_(17109.5)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'conductivity' dictionary
    if((dict.found("conductivity")) || (dict.dictName() == "conductivity"))
    {
        const dictionary& conductivityDict = dict.found("conductivity")?
        dict.subDict("conductivity"):
        dict;

        kappaInf_ = conductivityDict.lookupOrDefault<scalar>("kappaInf", 1.755);
        phi_ = conductivityDict.lookupOrDefault<scalar>("phi", 128.75e3);
        A0_ = conductivityDict.lookupOrDefault<scalar>("A0", 0.01926);
        Ax_ = conductivityDict.lookupOrDefault<scalar>("Ax", 1.06e-6);
        APu_ = conductivityDict.lookupOrDefault<scalar>("APu", 2.63e-8);
        B0_ = conductivityDict.lookupOrDefault<scalar>("B0", 2.39e-4);
        BPu_ = conductivityDict.lookupOrDefault<scalar>("BPu", 1.37e-13);
        D_ = conductivityDict.lookupOrDefault<scalar>("D", 5.27e9);
        E_ = conductivityDict.lookupOrDefault<scalar>("E", 17109.5);
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

Foam::conductivityUPuO2MagniMAMOX::~conductivityUPuO2MagniMAMOX()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUPuO2MagniMAMOX::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{
    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }
 
    const scalarField& Bu = Bu_->internalField();

    //- Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];              

        const scalar Bui = Bu[cellI];

        const scalar Ti = T[cellI];

        // Cut-off max porosity
        const scalar porI = min(porosity_[cellI], 0.95);

        // Number densities for heavy metals in fuel
        const scalar N_U = N_U_[cellI];
        const scalar N_Pu = N_Pu_[cellI];
        const scalar N_Am = (N_Am_ == nullptr) ? 0.0 : N_Am_->internalField()[cellI];
        const scalar N_Np = (N_Np_ == nullptr) ? 0.0 : N_Np_->internalField()[cellI];

        // Compute concentrations atoms/atomsHM
        const scalar conc_at_Pu = N_Pu / (N_Pu + N_Am + N_Np + N_U);    
        const scalar conc_at_Am = N_Am / (N_Pu + N_Am + N_Np + N_U);
        const scalar conc_at_Np = N_Np / (N_Pu + N_Am + N_Np + N_U);

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
        const scalar k_Np = 1.14;

        // Assign concentrations in mass/mass_fuel
        const scalar c_w_Pu = conc_at_Pu / k_Pu;
        const scalar c_w_Am = conc_at_Am / k_Am;
        const scalar c_w_Np = conc_at_Np / k_Np;

        // Compute deviation from stoichiometry
        const scalar x = 2 - OM_[cellI];

        // Computing Conductivity at fresh fuel k0            
        double term1 = 
            pow
            (
                A0_ + Ax_ * x + APu_ * c_w_Pu + AAm_ * c_w_Am + ANp_ * c_w_Np 
                + ( B0_ + BPu_ * c_w_Pu + BAm_ * c_w_Am + BNp_ * c_w_Np ) * Ti
                , 
                -1.0
            );
        double term2 = D_ / pow(Ti , 2.0) * exp(-E_ / Ti);
        double term3 = pow(1-porI , 2.5);
        double k0    = (term1 + term2) * term3;

        // Computing conductivity under irradiation (i.e. f(Bu))
        const scalar nominalValue = kappaInf_ + (k0 - kappaInf_)*exp(-Bui/phi_);

        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
