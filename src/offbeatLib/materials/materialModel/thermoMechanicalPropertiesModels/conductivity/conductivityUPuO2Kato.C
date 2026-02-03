
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

#include "conductivityUPuO2Kato.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(conductivityUPuO2Kato, 0);
    addToRunTimeSelectionTable
    (
        conductivityModel, 
        conductivityUPuO2Kato, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::conductivityUPuO2Kato::conductivityUPuO2Kato
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    conductivityModel(mesh, dict, baseDict, defaultModel),
    densityFrac_
    (
        dict.lookupOrDefault<scalar>
        (
            "densityFraction",
            baseDict.lookupOrDefault<scalar>("densityFraction", 0.945)
        )
    ),
    p_(nullptr),
    OM_(mesh_.lookupObject<volScalarField>("oxygenMetalRatio")),
    N_U_(mesh_.lookupObject<volScalarField>("N_U")),
    N_Pu_(mesh_.lookupObject<volScalarField>("N_Pu")),
    N_Am_(nullptr),
    N_Np_(nullptr),
    perturb(1)  
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'conductivity' dictionary
    if((dict.found("conductivity")) || (dict.dictName() == "conductivity"))
    {
        const dictionary& conductivityDict = dict.found("conductivity")?
        dict.subDict("conductivity"):
        dict;
        
        perturb = conductivityDict.lookupOrDefault<scalar>("perturb", 1.0);
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

Foam::conductivityUPuO2Kato::~conductivityUPuO2Kato()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::conductivityUPuO2Kato::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{

    scalar PiLim(1);
    scalar Pi(1);

    if (p_ == nullptr)
    {
        if ( mesh_.foundObject<volScalarField>("porosity") )
        {
            p_ = &mesh_.lookupObject<volScalarField>("porosity");
        }
    }

    //- Correct conductivity for every cell in addr
    forAll(addr, i)
    {   
        const label cellI = addr[i];

        // Temperature
        const scalar Ti = T[cellI];

        const scalar x = 2 - OM_[cellI];

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
        const scalar k_Am = 1.12;
        const scalar k_Np = 1.14;

        // Assign concentrations in mass/mass_fuel
        const scalar c_w_Am = conc_at_Am / k_Am;
        const scalar c_w_Np = conc_at_Np / k_Np;

        // Porosity
        if (p_ != nullptr) 
        {
            Pi = p_->internalField()[cellI];
        }
        
        // Maximum limit for porosity in this correlation, otherwise k->0 !
        if( !mesh_.foundObject<fvMesh>("referenceMesh") )
        {
            PiLim = min(Pi, 0.95); 
        }

        // Compute correction factor f(p) - Noascone
        const scalar fp = (1-PiLim)/(1+0.5*PiLim);

        // Compute denominator
        const scalar den = 1.595e-2 + 2.713*x + 3.583e-1*c_w_Am 
                         + 6.317e-2 *c_w_Np + (-2.625*x + 2.493)*1e-4*Ti;

        // Compute K(T)
        const scalar k = (1/(den)+1.541e11/pow(Ti,2.5)*exp(-1.522e4/Ti));

        // Correct k - Noascone method
        // const scalar nominalValue = fp * k;

        // Correct k - Barani method
        const scalar kappaHe(0.69);
        const scalar nominalValue = k * (kappaHe+2*k-2*Pi*(k-kappaHe))
                                      / (kappaHe+2*k + Pi*(k-kappaHe));
        
        // Assign k to scalar field
        sf[cellI] = nominalValue;
    }
}

// ************************************************************************* //
