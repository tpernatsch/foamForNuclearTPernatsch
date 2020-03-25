/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "diffusionNeutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionNeutronics, 0);

    addToRunTimeSelectionTable
    (
        neutronics,
        diffusionNeutronics,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionNeutronics::diffusionNeutronics
(
    fvMesh& mesh
)
:
    neutronics(mesh),//diffusionNeutronics is derived from neutronics
    xs_(mesh),
    Dalbedo_
    (
        IOobject
        (
            "Dalbedo",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimLength, 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    flux_(xs_.energyGroups()),
    fluxStar_(xs_.energyGroups()),
    prec_(xs_.precGroups()),
    precStar_(xs_.precGroups()),
    fluxStarAlbedo_
    (
        IOobject
        (
            "fluxStarAlbedo",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimArea/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    defaultFlux_
    (
        IOobject
        (
            "defaultFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    defaultPrec_
    (
        IOobject
        (
            "defaultPrec",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    oneGroupFlux_
    (
        IOobject
        (
            "oneGroupFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        defaultFlux_
    ),
    neutroSource_
    (
        IOobject
        (
            "neutroSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimVol/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    delayedNeutroSource_
    (
        IOobject
        (
            "delayedNeutroSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimVol/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    scatteringSourceExtra_
    (
        IOobject
        (
            "scatteringSourceExtra",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimless/dimVol/dimTime, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TFuel_
    (
        IOobject
        (
            "TFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TClad_
    (
        IOobject
        (
            "TClad",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TCool_
    (
        IOobject
        (
            "TCool",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoCool_
    (
        IOobject
        (
            "rhoCool",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, SMALL),
        zeroGradientFvPatchScalarField::typeName
    )
{
    #include "createNeutronicsFields.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusionNeutronics::~diffusionNeutronics()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusionNeutronics::correct
(
    scalar& residual, 
    label couplingIter
) 
{
    #include "solveNeutronics.H"
}

void Foam::diffusionNeutronics::getCouplingFieldRefs
(
    const objectRegistry& src,
    const meshToMesh& neutroToFluid
)
{
    //- Field names must reflect those defined in createCouplingFields.H
    TFuelOrig_ = 
        src.findObject<volScalarField>("bafflelessTFuelAv");
    TCladOrig_ = 
        src.findObject<volScalarField>("bafflelessTCladAv");
    TCoolOrig_ = 
        src.findObject<volScalarField>("bafflelessTCool");
    rhoCoolOrig_ = 
        src.findObject<volScalarField>("bafflelessRhoCool");
    if (liquidFuel_)
    {
        UOrig_ = 
            src.findObject<volVectorField>("bafflelessU");
        alphaOrig_ = 
            src.findObject<volScalarField>("bafflelessAlpha");
        alphatOrig_ = 
            src.findObject<volScalarField>("bafflelessAlphat");
        muOrig_ =
            src.findObject<volScalarField>("bafflelessMu");
    }
    else
    {
        UOrig_ = nullptr;
        alphaOrig_ = nullptr;
        alphatOrig_ = nullptr;
        muOrig_ = nullptr;
    }
}

void Foam::diffusionNeutronics::interpolateCouplingFields
(
    const meshToMesh& neutroToFluid
)
{
    neutroToFluid.mapTgtToSrc(*TFuelOrig_, plusEqOp<scalar>(), TFuel_);
    neutroToFluid.mapTgtToSrc(*TCladOrig_, plusEqOp<scalar>(), TClad_);
    neutroToFluid.mapTgtToSrc(*TCoolOrig_, plusEqOp<scalar>(), TCool_);
    neutroToFluid.mapTgtToSrc(*rhoCoolOrig_, plusEqOp<scalar>(), rhoCool_);
    if (liquidFuel_)
    {
        neutroToFluid.mapTgtToSrc(*UOrig_, plusEqOp<vector>(), UPtr_());
        neutroToFluid.mapTgtToSrc
        (
            *alphaOrig_, 
            plusEqOp<scalar>(), 
            alphaPtr_()
        );
        neutroToFluid.mapTgtToSrc
        (
            *alphatOrig_, 
            plusEqOp<scalar>(), 
            alphatPtr_()
        );
        neutroToFluid.mapTgtToSrc(*muOrig_, plusEqOp<scalar>(), muPtr_());
        phiPtr_() = fvc::flux(UPtr_());
        volScalarField diffCoeffOrig
        (
            (
                *alphatOrig_ 
            +   *muOrig_/xs_.ScNo()
            )/(*rhoCoolOrig_)
        ); 
        neutroToFluid.mapTgtToSrc
        (
            diffCoeffOrig, 
            plusEqOp<scalar>(), 
            diffCoeffPrecPtr_()
        );

        UPtr_().correctBoundaryConditions();
        alphaPtr_().correctBoundaryConditions();
        alphatPtr_().correctBoundaryConditions();
        diffCoeffPrecPtr_().correctBoundaryConditions();
    }

    TFuel_.correctBoundaryConditions();
    TClad_.correctBoundaryConditions();
    TCool_.correctBoundaryConditions();
    rhoCool_.correctBoundaryConditions();
}

// ************************************************************************* //
