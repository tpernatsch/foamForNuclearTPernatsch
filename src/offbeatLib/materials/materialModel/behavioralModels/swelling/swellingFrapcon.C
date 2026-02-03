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

#include "swellingFrapcon.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingFrapcon, 0);
    addToRunTimeSelectionTable(swellingModel, swellingFrapcon, dictionary);

    const char* swellingFrapcon::group_ = 
        ("behavioralModels::swelling::UO2Frapcon");

    defineParameter(swellingFrapcon, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingFrapcon, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingFrapcon::swellingFrapcon
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict  
)
:
    swellingModel(mesh, dict, baseDict),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    densityName_(dict.lookupOrDefault<word>("densityName", "rho" )),
    rho_(nullptr),
    par1(6000),
    par2(80000),
    par3(2.974e+10),
    par4(2.315e-23),
    par5(86.4),
    par6(3.211e-23)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'swelling' dictionary
    if((dict.found("swelling")) || (dict.dictName() == "swelling"))
    {
        const dictionary& swellingDict = dict.found("swelling")?
        dict.subDict("swelling"):
        dict;

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 6000);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", 80000);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", 2.974e+10);
        par4 = swellingDict.lookupOrDefault<scalar>("par4", 2.315e-23);
        par5 = swellingDict.lookupOrDefault<scalar>("par5", 86.4);
        par6 = swellingDict.lookupOrDefault<scalar>("par6", 3.211e-23);
    }  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingFrapcon::~swellingFrapcon()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingFrapcon::correct
(
    const scalarField& T, 
    const labelList& addr
)
{   
    // Check which F and delta should be used
    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() : swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() : swellingModel::delta_epsilonSwelling();

    const fvMesh& mesh_ = epsilonSwelling_.mesh();
    
    if(Bu_ == nullptr or rho_ == nullptr)
    {
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
        rho_ = &mesh_.lookupObject<volScalarField>(densityName_);
    }

    symmTensorField& swellingI = epsilonSwelling_.ref();
    
    const volScalarField& BuI = *Bu_;
    const volScalarField& rhoI = *rho_;
    
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
            
        const double a1 = par1;
        const double a2 = par2;
        const double b1 = par3 * par4 * par5;
        const double b2 = par3 * par6 * par5;
        
        const scalar Bui = BuI[cellI]/0.881;

        scalar nominalValue(0.0);
    
        if (Bui  < a1)
        {
                nominalValue = 0.0;
        }
        else if (Bui < a2)
        {
                nominalValue = (Bui - a1)*b1 * rhoI[cellI]/3;
        }
        else
        {
                nominalValue = ((a2 - a1)*b1 + (Bui - a2)*b2) * rhoI[cellI]/3;
        }
        
        // Perturb for sensitivity analysis
        nominalValue = nominalValue * F.value() + delta.value(); 

        swellingI[cellI] = nominalValue * I;
    }

    if(mesh_.foundObject<volSymmTensorField>("intragranularGasSwelling"))
    {
        const symmTensorField& Intragranular_gas_swelling = 
        mesh_.lookupObject<volSymmTensorField>("intragranularGasSwelling").internalField();

        const symmTensorField& Intergranular_gas_swelling = 
        mesh_.lookupObject<volSymmTensorField>("intergranularGasSwelling").internalField();

        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];
            swellingI[cellI] += Intragranular_gas_swelling[cellI]/3 
                              + Intergranular_gas_swelling[cellI]/3;
        }
    }

}


// ************************************************************************* //
