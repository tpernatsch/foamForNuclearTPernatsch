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

#include "swellingWrightShamHastelloyN.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingWrightShamHastelloyN, 0);
    addToRunTimeSelectionTable(swellingModel, swellingWrightShamHastelloyN, dictionary);

    const char* swellingWrightShamHastelloyN::group_ = 
        ("behavioralModels::swelling::WrightShamHastelloyN");

    defineParameter(swellingWrightShamHastelloyN, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingWrightShamHastelloyN, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingWrightShamHastelloyN::swellingWrightShamHastelloyN
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict
)
:
    swellingModel(mesh, dict, baseDict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    par1(9.845e-1),
    par2(4.385e-1),
    par3(-9.81e-1),
    par4(490.)
{  
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'swelling' dictionary
    if((dict.found("swelling")) || (dict.dictName() == "swelling"))
    {
        const dictionary& swellingDict = dict.found("swelling")?
        dict.subDict("swelling"):
        dict;

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 9.845e-1);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", 4.385e-1);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", -9.81e-1);
        par4 = swellingDict.lookupOrDefault<scalar>("par4", 490.);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingWrightShamHastelloyN::~swellingWrightShamHastelloyN()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingWrightShamHastelloyN::correct
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
    const symmTensorField& swellingIOld = epsilonSwelling_.oldTime().internalField();
    
    const scalarField& fastFluenceI = fastFluence_->internalField();
    const scalarField& fastFluenceIOld = fastFluence_->oldTime().internalField();

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        const scalar fastFluence = fastFluenceI[cellI]*1e-22*5;  //in dpa
        const scalar fastFluenceOld = fastFluenceIOld[cellI]*1e-22*5; //in dpa
        
        const scalar fluxfac    = par1*pow(fastFluence,par2) + par3;  //in dpa
        const scalar fluxfacOld = par1*pow(fastFluenceOld,par2) + par3; //in dpa
 
        
        const scalar Ti = T[cellI] - 273.15;
        const scalar sweexp = (Ti-par4)/100.0;
        const scalar tfact  = exp(-(sweexp*sweexp));
        
        const scalar nominalValue = tfact*fluxfac;
        const scalar nominalValueOld = tfact*fluxfacOld;
        
        
        scalar swellIncrement = 0.01/3.0*(nominalValue - nominalValueOld);
        
        // Perturb for sensitivity analysis
        swellIncrement = swellIncrement * F.value() + delta.value();

        swellingI[cellI] = swellingIOld[cellI];

        swellingI[cellI].xx() += (swellIncrement);
        swellingI[cellI].yy() += (swellIncrement);
        swellingI[cellI].zz() += (swellIncrement);
    }

}


// ************************************************************************* //
