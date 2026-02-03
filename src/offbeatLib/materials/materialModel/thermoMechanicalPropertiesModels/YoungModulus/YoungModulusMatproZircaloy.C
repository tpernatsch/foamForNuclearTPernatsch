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

#include "YoungModulusMatproZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(YoungModulusMatproZircaloy, 0);
    addToRunTimeSelectionTable
    (
        YoungModulusModel, 
        YoungModulusMatproZircaloy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::YoungModulusMatproZircaloy::YoungModulusMatproZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    YoungModulusModel(mesh, dict, baseDict, defaultModel),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence")),
    phi_(nullptr),
    cW_(dict.lookupOrDefault("coldWork", 0.0)),
    oxC_(dict.lookupOrDefault("oxyCon", 0.0)),
    par1(6.61e11),
    par2(5.912e8),
    par3(2.6e10),
    par4(0.88),
    par5(0.12),
    par6(1e25),
    par7(1.088e11),
    par8(5.475e7),
    par9(9.21e10),
    par10(4.05e7)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'YoungModulus' dictionary
    if((dict.found("YoungModulus")) || (dict.dictName() == "YoungModulus"))
    {
        const dictionary& YoungModulusDict = dict.found("YoungModulus")?
        dict.subDict("YoungModulus"):
        dict;

        par1 = YoungModulusDict.lookupOrDefault<scalar>("par1", 6.61e11);        
        par2 = YoungModulusDict.lookupOrDefault<scalar>("par2", 5.912e8);
        par3 = YoungModulusDict.lookupOrDefault<scalar>("par3", 2.6e10);
        par4 = YoungModulusDict.lookupOrDefault<scalar>("par4", 0.88);
        par5 = YoungModulusDict.lookupOrDefault<scalar>("par5", 0.12);
        par6 = YoungModulusDict.lookupOrDefault<scalar>("par6", 1e25);
        par7 = YoungModulusDict.lookupOrDefault<scalar>("par7", 1.088e11);
        par8 = YoungModulusDict.lookupOrDefault<scalar>("par8", 5.475e7);
        par9 = YoungModulusDict.lookupOrDefault<scalar>("par9", 9.21e10);
        par10 = YoungModulusDict.lookupOrDefault<scalar>("par10", 4.05e7);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::YoungModulusMatproZircaloy::~YoungModulusMatproZircaloy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::YoungModulusMatproZircaloy::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{ 
    if(phi_ == nullptr)
    {
        if(mesh_.foundObject<volScalarField>(fastFluenceName_))
        {
            phi_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
        }
    }

    forAll(addr, i)
    {   
        const label cellI = addr[i];        

        const scalar Ti = T[cellI];

        // Get neutron fluence and convert to n/m2
        scalar phi = 0.0;
        if(phi_ != nullptr)
        {
            phi  = phi_->internalField()[cellI] * 1e4;
        }

        // The validity range is changed by one degree on each side to allow 
        // for rounding errors
        if (Ti < 290 || Ti > 1801)
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", is out of range 290 < T < 1800 K";
        }

        scalar nominalValue(0.0);
        
        // Alpha phase
        if (Ti < 1073) 
        {
            // Oxygen concentration effect
            scalar K1 = (par1 + par2*Ti)*oxC_;
            
            // Cold work effect
            scalar K2 = -par3*cW_;
            
            // Fast neutron fluence effect 
            scalar K3 = par4 + par5*exp(-phi/par6);
            
            nominalValue = (par7 - par8*Ti + K1 + K2) / K3;
        }
        // Interpolate between alpha and beta values
        else if (Ti < 1273)
        {
            scalar TAlpha = 1073;
            scalar TBeta = 1273;
            
            // Oxygen concentration effect
            scalar K1 = (par1 + par2*TAlpha)*oxC_;
            // Cold work effect
            scalar K2 = -par3*cW_;       
            // Fast neutron fluence effect 
            scalar K3 = par4 + par5*exp(-phi/par6);
            
            scalar EAlpha = (par7 - par8*TAlpha + K1 + K2) / K3;
            scalar EBeta = par9 - par10*TBeta;
            
            nominalValue = EAlpha*(TBeta - Ti)/(TBeta - TAlpha) + EBeta*(Ti - TAlpha)/(TBeta - TAlpha);
        }
        // Beta phase
        else
        {
            nominalValue = (par9- par10*Ti);
        }

        sf[cellI] = nominalValue;
    }
}




// ************************************************************************* //
