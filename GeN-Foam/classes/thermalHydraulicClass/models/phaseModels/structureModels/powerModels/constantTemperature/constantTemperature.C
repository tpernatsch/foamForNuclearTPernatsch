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

#include "constantTemperature.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(constantTemperature, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        constantTemperature, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::constantTemperature::constantTemperature
(
    const structureModel& structure,
    const dictionary& dict,
    const word name
)
:
    powerModel
    (
        structure,
        dict,
        name
    ),
    T_
    (
        IOobject
        (
            "T."+typeName+"."+name,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("T", dimTemperature, dict),
        zeroGradientFvPatchScalarField::typeName
    )
{
    T_.primitiveFieldRef() *= cellField_;
    T_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::constantTemperature::~constantTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::powerModels::constantTemperature::T() const
{
    tmp<volScalarField> tT
    (
        new volScalarField(T_)
    );
    volScalarField& T = tT.ref();
    T.primitiveFieldRef() *= cellField_;
    T.correctBoundaryConditions();
    return tT;
}


// ************************************************************************* //
