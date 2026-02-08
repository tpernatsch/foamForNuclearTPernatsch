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

#include "PoissonRatioMatproZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PoissonRatioMatproZircaloy, 0);
    addToRunTimeSelectionTable
    (
        PoissonRatioModel, 
        PoissonRatioMatproZircaloy, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PoissonRatioMatproZircaloy::PoissonRatioMatproZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict,
    const word defaultModel
)
:
    PoissonRatioModel(mesh, dict, baseDict, defaultModel),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence")),
    YoungModulusName_(dict.lookupOrDefault<word>("YoungModulusName", "E")),
    phi_(nullptr),
    E_(nullptr),
    cW_(dict.lookupOrDefault("coldWork", 0.0)),
    oxC_(dict.lookupOrDefault("oxyCon", 0.0)),
    par1(7.07e11),
    par2(2.315e8),
    par3(2.6e10),
    par4(0.88),
    par5(0.12),
    par6(1e25),
    par7(4.04e10),
    par8(2.168e7),
    par9(3.49e10),
    par10(1.66e7)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'PoissonRatio' dictionary
    if((dict.found("PoissonRatio")) || (dict.dictName() == "PoissonRatio"))
    {
        const dictionary& PoissonRatioDict = dict.found("PoissonRatio")?
        dict.subDict("PoissonRatio"):
        dict;

        par1 = PoissonRatioDict.lookupOrDefault<scalar>("par1", 7.07e11);        
        par2 = PoissonRatioDict.lookupOrDefault<scalar>("par2", 2.315e8);
        par3 = PoissonRatioDict.lookupOrDefault<scalar>("par3", 2.6e10);
        par4 = PoissonRatioDict.lookupOrDefault<scalar>("par4", 0.88);
        par5 = PoissonRatioDict.lookupOrDefault<scalar>("par5", 0.12);
        par6 = PoissonRatioDict.lookupOrDefault<scalar>("par6", 1e25);
        par7 = PoissonRatioDict.lookupOrDefault<scalar>("par7", 4.04e10);
        par8 = PoissonRatioDict.lookupOrDefault<scalar>("par8", 2.168e7);
        par9 = PoissonRatioDict.lookupOrDefault<scalar>("par9", 3.49e10);
        par10 = PoissonRatioDict.lookupOrDefault<scalar>("par10", 1.66e7);
    }

}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PoissonRatioMatproZircaloy::~PoissonRatioMatproZircaloy()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::PoissonRatioMatproZircaloy::correct
(
    scalarField& sf, 
    const scalarField& T, 
    const labelList& addr
)
{     
    if(phi_ == nullptr or E_ == nullptr)
    {        
        if(mesh_.foundObject<volScalarField>(fastFluenceName_))
        {
            phi_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
        }

        E_ = &mesh_.lookupObject<volScalarField>(YoungModulusName_);
    }

    forAll(addr, i)
    {   
        const label cellI = addr[i];        

        const scalar Ti = T[cellI];
        
        // Get neutron fluence in n/cm2
        scalar phi = 0.0;
        if(phi_ != nullptr)
        {
            phi  = phi_->internalField()[cellI];
        }

        const scalar E = E_->internalField()[cellI];

        //- Shear modulus
        scalar G(0.0);

        //- Set the shear modulus, Matpro model
        //- The validity range is changed by one degree on each side to allow 
        //- for rounding errors
        if (Ti < 290 || Ti > 1801)
        {
            WarningInFunction
                << "Supplied temperature, " << Ti << ", is out of range 290 < T < 1800 K";
        }
       
        if (Ti < 1073) //- Alpha phase
        {
            // Oxygen concentration effect
            scalar K1 = (par1 - par2*Ti)*oxC_;
            
            // Cold work effect
            scalar K2 = -par3*cW_;
            
            // Fast neutron fluence effect 
            scalar K3 = par4 + par5*exp(-phi/par6);
            
            G = (par7 - par8*Ti + K1 + K2) / K3;
        }
        else if (Ti < 1273) //- Interpolate between alpha and beta values
        {
            scalar TAlpha = 1073;
            scalar TBeta = 1273;
            
            // Oxygen concentration effect
            scalar K1 = (par1 - par2*TAlpha)*oxC_;

            // Cold work effect
            scalar K2 = -par3*cW_;  

            // Fast neutron fluence effect 
            scalar K3 = par4 + par5*exp(-phi/par6);
            
            scalar GAlpha = (par7 - par8*TAlpha + K1 + K2) / K3;

            scalar GBeta = par9 - par10*TBeta;
            
            G = GAlpha*(TBeta - Ti)/(TBeta - TAlpha) + GBeta*(Ti - TAlpha)/(TBeta - TAlpha);
        }
        else //- Beta phase
        {
             G = (par9 - par10*Ti);
        }

        const scalar nominalValue = (E/(2*G) - 1);

        sf[cellI] = nominalValue;
    }
}

        
        


// ************************************************************************* //
