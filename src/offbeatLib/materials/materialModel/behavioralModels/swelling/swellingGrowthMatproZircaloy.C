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

#include "swellingGrowthMatproZircaloy.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingGrowthMatproZircaloy, 0);
    addToRunTimeSelectionTable(swellingModel, swellingGrowthMatproZircaloy, dictionary);

    const char* swellingGrowthMatproZircaloy::group_ = 
        ("behavioralModels::swelling::growthMatproZircaloy");

    defineParameter(swellingGrowthMatproZircaloy, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingGrowthMatproZircaloy, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingGrowthMatproZircaloy::swellingGrowthMatproZircaloy
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict 
)
:
    swellingModel(mesh, dict, baseDict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(mesh_.lookupObject<volScalarField>(fastFluenceName_)),
    cW_(dict.lookupOrDefault("coldWork", 0.0)),
    A(1.407e-16),
    par1(240.8),
    fz(0.05)
{  
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'swelling' dictionary
    if((dict.found("swelling")) || (dict.dictName() == "swelling"))
    {
        const dictionary& swellingDict = dict.found("swelling")?
        dict.subDict("swelling"):
        dict;

        A = swellingDict.lookupOrDefault<scalar>("A", 1.407e-16);
        par1 = swellingDict.lookupOrDefault<scalar>("par1", 240.8);
        fz = swellingDict.lookupOrDefault<scalar>("fz", 0.05);
    }  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingGrowthMatproZircaloy::~swellingGrowthMatproZircaloy()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingGrowthMatproZircaloy::correct
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

    symmTensorField& swellingI = epsilonSwelling_.ref();

    const scalarField& oldTimeSwellI = 
    epsilonSwelling_.oldTime().internalField();

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        const scalar Ti = T[cellI];

        //- Scale fast fluence from n/cm2/s to n/m2/s
        const scalar fastFluencei = fastFluence_[cellI]*1e4;


        scalar temp(0.0);
                
        if(Ti < 40 + 273.15)
        {
            temp = 40 + 273.15;
        }
        else if(Ti < 360 + 273.15)
        {
            temp = Ti;
        }
        else
        {
            temp = 360 + 273.15;
        }
        
        const scalar nominalValue = 
        A*exp(par1/temp)*pow(fastFluencei, 0.5)*(1 - 3*fz)*(1 + 2*cW_);  

        // Perturb for sensitivity analysis
        nominalValue = nominalValue * F.value() + delta.value(); 

        swellingI[cellI] = oldTimeSwellI[cellI];  

        swellingI[cellI] += nominalValue;
    }
}


// ************************************************************************* //
