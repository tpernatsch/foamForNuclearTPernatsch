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

#include "swellingFeCrAl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingFeCrAl, 0);
    addToRunTimeSelectionTable(swellingModel, swellingFeCrAl, dictionary);
    
    const char* swellingFeCrAl::group_ = 
        ("behavioralModels::swelling::FeCrAl");

    defineParameter(swellingFeCrAl, F_epsilonSwelling_, 
        "F_epsilonSwelling", (dimless), 1.0);
        
    defineParameter(swellingFeCrAl, delta_epsilonSwelling_, 
        "delta_epsilonSwelling", (dimless), 0.0);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingFeCrAl::swellingFeCrAl
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& baseDict   
)
:
    swellingModel(mesh, dict, baseDict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    par1(4.5E-29)
{
    // This part is kept only for retro-compatibility.
    // In the future, parameters are either eliminated or looked for by default 
    // in the 'swelling' dictionary
    if((dict.found("swelling")) || (dict.dictName() == "swelling"))
    {
        const dictionary& swellingDict = dict.found("swelling")?
        dict.subDict("swelling"):
        dict;

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 4.5E-29);
    }  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingFeCrAl::~swellingFeCrAl()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingFeCrAl::correct
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
        const fvMesh& mesh_ = epsilonSwelling_.mesh();
        
        fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
    }

    symmTensorField& swellingI = epsilonSwelling_.ref();
    
    const volScalarField& phii = *fastFluence_;
    
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        //- The model requires n/m^2
        const scalar phi = phii[cellI]*1e4;

        scalar nominalValue(par1*phi);

        // Perturb for sensitivity analysis
        nominalValue = nominalValue * F.value() + delta.value(); 

        swellingI[cellI] = nominalValue*I;
    }
}


// ************************************************************************* //
