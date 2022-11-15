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

#include "neutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neutronics, 0);
    defineRunTimeSelectionTable(neutronics, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neutronics::neutronics
(
    fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "neutronicsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    reactorState_
    (
        IOobject
        (
            "reactorState",
            mesh_.time().timeName(),
            "uniform",
            mesh_.time(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        )
    ),
    keff_(reactorState_.lookupOrDefault("keff",1.0)),
    pTarget_(reactorState_.lookupOrDefault("pTarget",1.0)),
    powerDensity_
    (
        IOobject
        (
            "powerDensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    secondaryPowerDenisty_
    (
        IOobject
        (
            "secondaryPowerDenisty",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    disp_
    (
        IOobject
        (
            "disp",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("d_zero", dimLength, vector::zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    initialResidual_(1.0),
    eigenvalueNeutronics_
    (
        IOdictionary::lookupOrDefault("eigenvalueNeutronics", false)
    ),
    liquidFuel_
    (
        mesh.time().controlDict().lookupOrDefault("liquidFuel", false)
    )
{
    Info << "Initial keff = " << keff_ << endl;
}

// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::neutronics> Foam::neutronics::New
(
    fvMesh& mesh
)
{
    word modelName;

    IOdictionary dict
    (
        IOobject
        (
            "neutronicsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dict.lookup("model") >> modelName;

    Info<< "Selecting neutronics model type " << modelName << endl;

    auto* ctorPtr = dictionaryConstructorTable(modelName);

    if (!ctorPtr)
    {
        FatalErrorIn
        (
            "neutronics::New(const volScalarField&, "
            "const volVectorField&, basicThermo&)"
        )   << "Unknown neutronics model " << modelName
            << endl << endl
            << "Valid models types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }
    return
        autoPtr<neutronics>
        (
            ctorPtr(mesh)
        );
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::neutronics::deformMesh
(
    const meshToMesh& TMToNeutro,
    const volVectorField& dispOrig
)
{
    const volPointInterpolation& neutroMeshPointInterpolation =
        volPointInterpolation::New(mesh_);

    tmp<pointVectorField> neutroPointsDisplacementOld =
        neutroMeshPointInterpolation.interpolate(disp_);

    disp_ *= 0.0;
    TMToNeutro.mapSrcToTgt(dispOrig, plusEqOp<vector>(), disp_);
    disp_.correctBoundaryConditions();

    tmp<pointVectorField> neutroPointsDisplacement =
        neutroMeshPointInterpolation.interpolate(disp_);

    tmp<pointField> displacedPoints =
        mesh_.points()
    +   neutroPointsDisplacement->internalField()
    -   neutroPointsDisplacementOld->internalField();

    mesh_.movePoints(displacedPoints);
}

// ************************************************************************* //
