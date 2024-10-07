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

#include "swellingMATPRO.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(swellingMATPRO, 0);
    addToRunTimeSelectionTable(swellingModel, swellingMATPRO, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swellingMATPRO::swellingMATPRO
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    swellingModel(mesh, dict),
    burnupName_(dict.lookupOrDefault<word>("burnupName", "Bu" )),
    Bu_(nullptr),
    densityName_(dict.lookupOrDefault<word>("densityName", "rho" )),
    rho_(nullptr),
    par1(5.577e-5),
    par2(1.96e-31),
    par3(2800),
    par4(11.73),
    par5(-0.0162),
    par6(-0.0178),
    perturb(1.0)
{
    if(dict.found("swelling"))
    {
        const dictionary& swellingDict = dict.subDict("swelling");

        par1 = swellingDict.lookupOrDefault<scalar>("par1", 5.577e-5);
        par2 = swellingDict.lookupOrDefault<scalar>("par2", 1.96e-31);
        par3 = swellingDict.lookupOrDefault<scalar>("par3", 2800);
        par4 = swellingDict.lookupOrDefault<scalar>("par4", 11.73);
        par5 = swellingDict.lookupOrDefault<scalar>("par5", -0.0162);
        par6 = swellingDict.lookupOrDefault<scalar>("par6", -0.0178);

        perturb = swellingDict.lookupOrDefault<scalar>("perturb", 1.0);
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::swellingMATPRO::~swellingMATPRO()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::swellingMATPRO::correct
(
    const scalarField& T,
    const labelList& addr
)
{
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

        swellingI[cellI] = swellingIOld[cellI]
                          +(solidSwellIncrement+gasSwellIncrement)/3.0*perturb*I;
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
