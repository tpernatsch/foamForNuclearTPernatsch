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

#include "FFPair.H"
#include "FFDragFactor.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FFDragFactor, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FFDragFactor::FFDragFactor
(
    const FFPair& pair,
    const dictionary& dict
)
:
    IOdictionary
    (
        IOobject
        (
            "FFDragFactor."+typeName+"."+pair.name()+"."+dict.dictName(),
            pair.mesh().time().timeName(),
            pair.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    pair_(pair),
    coeffReg_(*this),
    rhoc_(pair_.rhoContinuous()),
    Dhd_(pair_.DhDispersed()),
    magUr_(pair_.magUr())
{
    fdPtr_.reset
    (
        FFDragCoefficientModel::New
        (
            pair_,
            *this,
            coeffReg_
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FFDragFactor::correctField(volScalarField& Kd) const
{
    const volScalarField& a1(pair_.fluid1());
    const volScalarField& a2(pair_.fluid2());
    forAll(pair_.mesh().cells(), i)
    {
        Kd[i] = 
            (a1[i]+a2[i])*(0.5/Dhd_[i])*rhoc_[i]*magUr_[i]*fdPtr_->value(i);
    }
}


// ************************************************************************* //
