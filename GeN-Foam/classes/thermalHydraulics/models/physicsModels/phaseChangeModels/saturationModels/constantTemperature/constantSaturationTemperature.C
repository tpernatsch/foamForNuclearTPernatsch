/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "constantSaturationTemperature.H"
#include "fluid.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{
    defineTypeNameAndDebug(constantTemperature, 0);
    addToRunTimeSelectionTable
    (
        saturationModel,
        constantTemperature,
        saturationModels
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::saturationModels::constantTemperature::
constantTemperature
(
    const dictionary& dict, 
    const fluid& fluid1,
    const fluid& fluid2
)
:
    saturationModel
    (
        dict,
        fluid1,
        fluid2
    ),
    Tsat_("T", dimTemperature, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantTemperature::pSat
(
    const volScalarField& T
) const
{
    return fluid1_.thermo().p();
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantTemperature::pSatPrime
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "pSatPrime",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimPressure/dimTemperature, 0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantTemperature::lnPSat
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            log(fluid1_.thermo().p())
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::saturationModels::constantTemperature::Tsat
(
    const volScalarField& p
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Tsat",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            Tsat_
        )
    );
}


// ************************************************************************* //
