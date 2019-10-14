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

#include "constantFluidPower.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace powerModels
{
    defineTypeNameAndDebug(constantFluidPower, 0);
    addToRunTimeSelectionTable
    (
        powerModel, 
        constantFluidPower, 
        powerModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::powerModels::constantFluidPower::constantFluidPower
(
    const structureModel& structure,
    const dictionary& dict,
    const word name,
    volScalarFieldTable& fields
)
:
    powerModel
    (
        structure,
        dict,
        name,
        fields
    ),
    P_
    (
        this->initInFieldsTable
        (
            volScalarField
            (
                IOobject
                (
                    typeName+".P",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("P", dimPower/dimVol, dict),
                zeroGradientFvPatchScalarField::typeName
            )
        )
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::powerModels::constantFluidPower::~constantFluidPower()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> 
Foam::powerModels::constantFluidPower::heatSource
(
    const volScalarField& Tfluid,
    const volScalarField& h
) const
{
    return 
        tmp<volScalarField>
        (
            new volScalarField
            (
                P_
            )
        );
}


// ************************************************************************* //
