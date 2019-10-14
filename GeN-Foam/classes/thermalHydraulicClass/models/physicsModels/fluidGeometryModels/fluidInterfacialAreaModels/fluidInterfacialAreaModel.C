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

#include "fluidInterfacialAreaModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fluidInterfacialAreaModel, 0);
    defineRunTimeSelectionTable
    (
        fluidInterfacialAreaModel, 
        fluidInterfacialAreaModels
    );
}

const Foam::dimensionSet 
    Foam::fluidInterfacialAreaModel::dimIA(0, -1, 0, 0, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModel::fluidInterfacialAreaModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& fluid1,
    const fluid& fluid2
)
:
    IOdictionary
    (
        IOobject
        (
            typeName,
            fluid1.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    fluid1_(fluid1),
    fluid2_(fluid2),
    residualAlpha_
    (
        Foam::sqrt
        (
            (
                fluid1_.residualAlpha()
            +   fluid2_.residualAlpha()
            )/2
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fluidInterfacialAreaModel::~fluidInterfacialAreaModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
