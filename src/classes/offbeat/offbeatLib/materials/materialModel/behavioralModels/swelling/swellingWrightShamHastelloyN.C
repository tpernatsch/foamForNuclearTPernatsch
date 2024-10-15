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
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingWrightShamHastelloyN::swellingWrightShamHastelloyN
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    swellingModel(mesh, dict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    par1(9.845e-1),
    par2(4.385e-1),
    par3(-9.81e-1),
    par4(490.),
    perturb(1.0)
{
    if(dict.found("swelling"))
    {
        const dictionary& swellingDict = dict.subDict("swelling");

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 9.845e-1);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", 4.385e-1);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", -9.81e-1);
        par4 = swellingDict.lookupOrDefault<scalar>("par4", 490.);


        perturb = swellingDict.lookupOrDefault<scalar>("perturb", 1.0);
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
        
        
        const scalar swellIncrement = 0.01/3.0*(nominalValue - nominalValueOld);
        
        swellingI[cellI] = swellingIOld[cellI];

        swellingI[cellI].xx() += (swellIncrement*perturb);
        swellingI[cellI].yy() += (swellIncrement*perturb);
        swellingI[cellI].zz() += (swellIncrement*perturb);
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
