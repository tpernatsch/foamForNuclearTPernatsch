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

#include "UPuO2meltingMagni2020.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(UPuO2meltingMagni2020, 0);
    addToRunTimeSelectionTable(failureModel, UPuO2meltingMagni2020, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::UPuO2meltingMagni2020::UPuO2meltingMagni2020
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    failureModel(mesh, dict, baseDict),
    burnupName_(failureModelDict_.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    TmUO2_(failureModelDict_.lookupOrDefault<scalar>("Tmelt", 3147.0)),
    OM_(mesh_.lookupObject<volScalarField>("oxygenMetalRatio")),
    N_U_(mesh_.lookupObject<volScalarField>("N_U")),
    N_Pu_(mesh_.lookupObject<volScalarField>("N_Pu")),
    N_Am_(mesh_.lookupObject<volScalarField>("N_Am")),
    N_Np_(mesh_.lookupObject<volScalarField>("N_Np")),
    gammaPu_(364.85),
    gammaX_(1014.15),
    TmInf_(2964.92),
    delta_(40.43e3)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UPuO2meltingMagni2020::~UPuO2meltingMagni2020()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::UPuO2meltingMagni2020::checkFailure
(
    const labelList& addr
)
{
    if(Bu_ == nullptr)
    {        
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
    }
 
    const scalarField& Bu = Bu_->internalField();

    // Set failure switch to false
    failed_ = false;

    // Set counter of failed cells to 0
    scalar counter(0);

    
    // Constant reference to the T field
    const volScalarField& T(mesh_.lookupObject<volScalarField>("T"));

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        // Number densities for heavy metals in fuel
        const scalar N_U = N_U_[cellI];
        const scalar N_Pu = N_Pu_[cellI];
        const scalar N_Am = N_Am_[cellI];
        const scalar N_Np = N_Np_[cellI];

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

        // Compute deviation from stoichiometry
        const scalar x = 2 - OM_[cellI];

        // Compute Tm,0
        const scalar Tm0 = TmUO2_ - gammaPu_*c_w_Pu - gammaX_*x;

        // Compute Tm,irr
        const scalar TmIrr 
            = TmInf_ + (Tm0 - TmInf_)*Foam::exp(-Bu[cellI]/delta_);

        if ( T[cellI] >= TmIrr)
        {
            if( !mesh_.foundObject<volScalarField>("failedMaterial"))
            {
                //- Initialize failedMaterial field (from parent class "failureModel") 
                initializeFailedMaterialField();
            }

            failedMaterial_()[cellI] = 1;
            failed_ = true;
            counter++;
        }
    }

    if ( failed_ )
    {

        if ( stopIfFailed_ )
        {
            // Write fields files in the time-step directory at which the failure occurs
            mesh_.time().write();

            // Fatal error
            FatalErrorIn("Foam::UPuO2meltingMagni2020::checkFailure(const labelList& addr)") << nl
            << "    Failure occurred in " << counter << " out of " << addr.size() 
            << " cells of the material."  << nl      << exit(FatalError); 
        }

        //- Warning
        WarningIn("Foam::UPuO2meltingMagni2020::checkFailure(const labelList& addr)") << nl  
        << "    Failure occurred in " << counter << " out of " << addr.size() 
        << " cells of the material."  << nl      << endl;
    }
}


// ************************************************************************* //

