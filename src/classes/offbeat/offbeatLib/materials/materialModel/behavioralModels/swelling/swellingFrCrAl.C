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

#include "swellingFrCrAl.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingFrCrAl, 0);
    addToRunTimeSelectionTable(swellingModel, swellingFrCrAl, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingFrCrAl::swellingFrCrAl
(
    const fvMesh& mesh,
    const dictionary& dict   
)
:
    swellingModel(mesh, dict),
    fastFluenceName_(dict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluence_(nullptr),
    par1(4.5E-29),
    perturb(1.0)
{
    if(dict.found("swelling"))
    {
        const dictionary& swellingDict = dict.subDict("swelling");

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 4.5E-29);

        perturb = swellingDict.lookupOrDefault<scalar>("perturb", 1.0);
    }  
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingFrCrAl::~swellingFrCrAl()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingFrCrAl::correct
(
    const scalarField& T, 
    const labelList& addr
)
{     
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

        swellingI[cellI] = (nominalValue*perturb)*I;
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
