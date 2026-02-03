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

#include "swellingGrowthAIM11515Ti.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingGrowthAIM11515Ti, 0);
    addToRunTimeSelectionTable(swellingModel, swellingGrowthAIM11515Ti, dictionary);

    const char* swellingGrowthAIM11515Ti::group_ = 
        ("behavioralModels::swelling::growthAIM11515Ti");

    defineParameter(swellingGrowthAIM11515Ti, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingGrowthAIM11515Ti, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingGrowthAIM11515Ti::swellingGrowthAIM11515Ti
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    swellingModel(mesh, dict, baseDict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    par1(1.3e-5),
    par2(-2.5),
    par3(490),
    par4(100),
    par5(3.9)
{  
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'swelling' dictionary
    if((dict.found("swelling")) || (dict.dictName() == "swelling"))
    {
        const dictionary& swellingDict = dict.found("swelling")?
        dict.subDict("swelling"):
        dict;

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 1.3e-5);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", -2.5);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", 490);
        par4 = swellingDict.lookupOrDefault<scalar>("par4", 100);
        par5 = swellingDict.lookupOrDefault<scalar>("par5", 3.9);
    }  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingGrowthAIM11515Ti::~swellingGrowthAIM11515Ti()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingGrowthAIM11515Ti::correct
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
    
    if(fastFluence_ == nullptr)
    {        
        fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
    }

    symmTensorField& swellingI = epsilonSwelling_.ref();
    const symmTensorField& swellingIOld = 
        epsilonSwelling_.oldTime().internalField();

    const scalarField& fastFluenceI = fastFluence_->internalField();
    const scalarField& fastFluenceOldI = fastFluence_->oldTime().internalField();

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        const scalar Ti = T[cellI] - 273.15;

        //- Divide fast fluence by 1e22 n cm^-2
        const scalar fastFluence = fastFluenceI[cellI]*1e-22;
        const scalar fastFluenceOld = fastFluenceOldI[cellI]*1e-22;
        
        scalar nominalValue = 1.0/3.0 *
        par1 * exp(par2*(pow((Ti-par3)/par4,2)))*pow(fastFluence,par5)*1e-2;

        const scalar nominalValueOld = 1.0/3.0 *
        par1 * exp(par2*(pow((Ti-par3)/par4,2)))*pow(fastFluenceOld,par5)*1e-2;

        // Compute swelling increment wrt old time
        scalar swellIncrement(nominalValue-nominalValueOld);

        // Perturb for sensitivity analysis
        swellIncrement = swellIncrement * F.value() + delta.value(); 

        swellingI[cellI] = swellingIOld[cellI];

        swellingI[cellI].xx() += (swellIncrement);
        swellingI[cellI].yy() += (swellIncrement);
        swellingI[cellI].zz() += (swellIncrement);
    }
}


// ************************************************************************* //
