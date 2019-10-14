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

#include "structureInterfacialAreaPartitionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(structureInterfacialAreaPartitionModel, 0);
    defineRunTimeSelectionTable
    (
        structureInterfacialAreaPartitionModel, 
        structureInterfacialAreaPartitionModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::structureInterfacialAreaPartitionModel::
structureInterfacialAreaPartitionModel
(
    const objectRegistry& objReg,
    const dictionary& dict,
    const fluid& fluid1,
    const fluid& fluid2,
    const structureModel& structure
)
:
    IOdictionary
    (
        IOobject
        (
            "structureInterfacialAreaPartitionModel",
            fluid1.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    mesh_(fluid1.mesh()),
    fluid1_(fluid1),
    fluid2_(fluid2),
    structure_(structure)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::structureInterfacialAreaPartitionModel::~structureInterfacialAreaPartitionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// ************************************************************************* //
