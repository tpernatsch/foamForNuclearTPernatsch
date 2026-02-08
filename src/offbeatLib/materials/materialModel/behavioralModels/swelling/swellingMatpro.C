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

#include "swellingMatpro.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingMatpro, 0);
    addToRunTimeSelectionTable(swellingModel, swellingMatpro, dictionary);

    const char* swellingMatpro::group_ = 
        ("behavioralModels::swelling::UO2Matpro");

    defineParameter(swellingMatpro, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingMatpro, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingMatpro::swellingMatpro
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
    par1(5.577e-5),
    par2(1.96e-31),
    par3(2800),
    par4(11.73),
    par5(-0.0162),
    par6(-0.0178)
{  
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'swelling' dictionary
    if((dict.found("swelling")) || (dict.dictName() == "swelling"))
    {
        const dictionary& swellingDict = dict.found("swelling")?
        dict.subDict("swelling"):
        dict;

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 5.577e-5);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", 1.96e-31);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", 2800);
        par4 = swellingDict.lookupOrDefault<scalar>("par4", 11.73);
        par5 = swellingDict.lookupOrDefault<scalar>("par5", -0.0162);
        par6 = swellingDict.lookupOrDefault<scalar>("par6", -0.0178);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingMatpro::~swellingMatpro()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingMatpro::correct
(
    const scalarField& T,
    const labelList& addr
)
{
   // Check which F and delta should be used
    const dimensionedScalar& F =
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() :
    swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta =
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() :
    swellingModel::delta_epsilonSwelling();

    const fvMesh& mesh_ = epsilonSwelling_.mesh();

    if(Bu_ == nullptr or rho_ == nullptr)
    {
        Bu_ = &mesh_.lookupObject<volScalarField>(burnupName_);
        rho_ = &mesh_.lookupObject<volScalarField>(densityName_);
    }

    const symmTensorField& swellingIOld = epsilonSwelling_.oldTime().internalField();
    symmTensorField& swellingI = epsilonSwelling_.ref();

    //
    const scalarField& BuI = Bu_ -> internalField();
    const scalarField& BuIOld = Bu_ -> oldTime().internalField();
    const scalarField& rhoI = rho_ -> internalField();

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        // Convert Bu in MWd/tUO2 to FIMA
        const scalar Bu = BuI[cellI]*1e-3/0.881/937.06;
        const scalar BuOld = BuIOld[cellI]*1e-3/0.881/937.06;
        const scalar rho = rhoI[cellI];

        const scalar solidSwellIncrement = par1*rho*(Bu-BuOld);
        const scalar gasSwellIncrement = par2*rho*(Bu-BuOld)*pow(par3-T[cellI],par4)
                                        *exp(par5*(par3-T[cellI]))*exp(par6*rho*Bu);

        scalar nominalValue = (solidSwellIncrement+gasSwellIncrement)/3.0;

        // Perturb for sensitivity analysis
        nominalValue = nominalValue * F.value() + delta.value();

        swellingI[cellI] = swellingIOld[cellI] + nominalValue * I;
    }
}


// ************************************************************************* //
