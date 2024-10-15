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

#include "swellingGrowthGeneralized1515Ti.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingGrowthGeneralized1515Ti, 0);
    addToRunTimeSelectionTable(swellingModel, swellingGrowthGeneralized1515Ti, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingGrowthGeneralized1515Ti::swellingGrowthGeneralized1515Ti
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    swellingModel(mesh, dict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    par1(1.5e-3),
    par2(-2.5),
    par3(450),
    par4(100),
    par5(2.75),
    perturb(1.0)
{    
    if(dict.found("swelling"))
    {
        const dictionary& swellingDict = dict.subDict("swelling");

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 1.5e-3);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", -2.5);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", 450);
        par4 = swellingDict.lookupOrDefault<scalar>("par4", 100);
        par5 = swellingDict.lookupOrDefault<scalar>("par5", 2.75);

        perturb = swellingDict.lookupOrDefault<scalar>("perturb", 1.0);
    }  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingGrowthGeneralized1515Ti::~swellingGrowthGeneralized1515Ti()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingGrowthGeneralized1515Ti::correct
(
    const scalarField& T, 
    const labelList& addr
)
{
    if(fastFluence_ == nullptr)
    {        
        fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
    }

    symmTensorField& swellingI = epsilonSwelling_.ref();

    const scalarField& fastFluenceI = fastFluence_->internalField();

    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        
        const scalar Ti = T[cellI] - 273.15;

        //- Divide fast fluence by 1e22 n cm^-2
        const scalar fastFluence = fastFluenceI[cellI]*1e-22;
        
        const scalar nominalValue = 1.0/3.0 *
        par1 * exp(par2*(pow((Ti-par3)/par4,2)))*pow(fastFluence,par5)*1e-2;

        swellingI[cellI] = nominalValue*perturb*I;
    }

    // Check which F and delta should be used
    const dimensionedScalar& F = 
    isInUserParameters(this->group(), "F_epsilonSwelling") ?
    this->F_epsilonSwelling() : 
    swellingModel::F_epsilonSwelling();

    const dimensionedScalar& delta = 
    isInUserParameters(this->group(), "delta_epsilonSwelling") ?
    this->delta_epsilonSwelling() : 
    swellingModel::delta_epsilonSwelling();

    // Perturb epsilon swelling for sensitivity analysis
    applyParameters(epsilonSwelling_.ref(), addr, F, delta);

}


// ************************************************************************* //
