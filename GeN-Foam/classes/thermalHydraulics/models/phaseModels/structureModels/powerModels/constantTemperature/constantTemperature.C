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
#include "structureModel.H"
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
    const dictionary& dicts
)
:
    powerModel
    (
        structure,
        dicts
    ),
    T_
    (
        IOobject
        (
            "T."+typeName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE //AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("T", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
{
    structure_.setRegionField(this, T_, "T");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::constantTemperature::~constantTemperature()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::powerModels::constantTemperature::correctT(volScalarField& T) const
{
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        T[celli] = T_[celli];
    }
}

void Foam::powerModels::constantTemperature::powerOff()
{
    //- If you set iA to 0, the energy contribution from this powerModel to the
    //  fluid energy equation will be 0, equivalent to a "power" off scenario
    iA_ *= 0.0;
}

// ************************************************************************* //
