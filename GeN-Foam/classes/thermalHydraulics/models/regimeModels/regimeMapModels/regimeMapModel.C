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

#include "regimeMapModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regimeMapModel, 0);
    defineRunTimeSelectionTable
    (
        regimeMapModel, 
        regimeMapModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regimeMapModel::regimeMapModel
(
    const fvMesh& mesh,
    const dictionary& dict,
    const dictionary& physicsModelsDict
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName("regimeMapModel", typeName),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),    
    mesh_(mesh),
    physicsModelsDict_(physicsModelsDict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regimeMapModel::setRequiresModelCorrection()
{
    forAllIter
    (
        regimeTable,
        regimes_,
        iter
    )
    {
        regime& regime(iter()());
        if (!regime.isInterpolated())
        {
            regime.requiresModelCorrection() = regime.isCurrentlyPresent();
        }
        else
        {
            regime.regime1().requiresModelCorrection() = true;
            regime.regime2().requiresModelCorrection() = true;
        }
    }
}

// ************************************************************************* //
