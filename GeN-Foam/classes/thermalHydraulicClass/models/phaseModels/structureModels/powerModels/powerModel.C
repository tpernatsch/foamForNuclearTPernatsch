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

#include "powerModel.H"
#include "structureModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(powerModel, 0);
    defineRunTimeSelectionTable
    (
        powerModel, 
        powerModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModel::powerModel
(
    const structureModel& structure,
    const dictionary& dict,
    const word name
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(typeName, name),
            structure.mesh().time().timeName(),
            structure.mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    structure_(structure),
    mesh_(structure.mesh()),
    name_(name),
    cellList_(structure.cellLists()[name]),
    cellField_(structure.cellFields()[name]),
    iA_
    (
        IOobject
        (
            IOobject::groupName("iA", name),
            mesh_.time().timeName(),
            mesh_
        ),
        mesh_,
        dimensionedScalar::lookupOrDefault
        (
            "iA",
            dict,
            dimArea/dimVol,
            0
        ),
        zeroGradientFvPatchScalarField::typeName
    )
{    
    iA_.primitiveFieldRef() *= cellField_;
    iA_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
