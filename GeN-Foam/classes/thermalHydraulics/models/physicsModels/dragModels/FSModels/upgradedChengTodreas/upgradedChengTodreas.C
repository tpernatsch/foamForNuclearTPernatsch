/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "upgradedChengTodreas.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace dragModels
{
    defineTypeNameAndDebug(upgradedChengTodreas, 0);
    addToRunTimeSelectionTable(dragModel, upgradedChengTodreas, FSDragModels);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dragModels::upgradedChengTodreas::upgradedChengTodreas
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const FSPair& FSPair
)
:
    dragModel
    (
        objReg,
        dict,
        FSPair
    ),
    PbyD_(0),
    aL_(this->get<scalar>("aL")),
    bL_(this->get<scalar>("bL")),
    cL_(this->get<scalar>("cL")),
    aT_(this->get<scalar>("aT")),
    bT_(this->get<scalar>("bT")),
    cT_(this->get<scalar>("cT")),
    ReL_(0),
    ReT_(0)
{
    scalar P(this->get<scalar>("pinPitch"));
    scalar D(this->get<scalar>("pinDiameter"));
    PbyD_ = P/D;
    scalar PDm1(PbyD_-1.0);
    ReL_ = 320*(Foam::pow(10, PDm1));
    ReT_ = 10000*(Foam::pow(10, 0.7*(PDm1)));
    logReTReL_ = Foam::log(ReT_/ReL_);
    CL_ = aL_ + bL_*PDm1+cL_*sqr(PDm1);
    CT_ = aT_ + bT_*PDm1+cT_*sqr(PDm1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::dragModels::upgradedChengTodreas::fd(const scalar& Re) const
{
    if (Re < ReL_)
    {
        return CL_/Re;
    }
    else if (Re > ReT_)
    {
        return CT_/pow(Re, 0.18);
    }
    else 
    {
        scalar psi(Foam::log(Re/ReL_)/logReTReL_);
        return 
            pow(1.0-psi, 0.3333)*(1.0-pow(psi, 7))*CL_/Re*
        +   pow(psi, 0.3333)*CT_/pow(Re, 0.18);
    }
}

void Foam::dragModels::upgradedChengTodreas::correctKd(volTensorField& Kd) const
{    
    #include "calcKdFromFd.H"
}


// ************************************************************************* //
