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
    xs(mesh),
    flux(xs.energyGroups()),
    fluxStar(xs.energyGroups()),
    prec(xs.precGroups()),
    fluxStarAlbedo
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
        dimensionedScalar("", dimensionSet(0,-2,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    defaultFlux
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
    defaultPrec
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
        dimensionedScalar("", dimensionSet(0,-3,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Dalbedo
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
        dimensionedScalar("", dimensionSet(0,1,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    oneGroupFlux
    (
        IOobject
        (
            "oneGroupFlux",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        defaultFlux
    ),
    neutroSource
    (
        IOobject
        (
            "neutroSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    delayedNeutroSource
    (
        IOobject
        (
            "delayedNeutroSource",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    scatteringSourceExtra
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
        dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    U
    (
        IOobject
        (
            "U",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("", dimensionSet(0,1,-1,0,0,0,0), vector(0.0,0.0,0.0)),
        zeroGradientFvPatchVectorField::typeName
    ),
    phi
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U)
    ),
    porosity
    (
        IOobject
        (
            "porosity",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0),1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tfuel
    (
        IOobject
        (
            "Tfuel",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Tclad
    (
        IOobject
        (
            "Tclad",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoCool
    (
        IOobject
        (
            "rhoCool",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(1,-3,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TCool
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
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    diffCoeffPrec
    (
        IOobject
        (
            "diffCoeffPrec",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,2,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    )

{
    #include "createNeutronicsFields.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::diffusionNeutronics::~diffusionNeutronics()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::diffusionNeutronics::correct(const label couplingIter, scalar& residual, const bool& liquidFuel ) 
{
    #include "solveNeutronics.H"
}


void Foam::diffusionNeutronics::getFields(
    const volScalarField& TfuelOrig, 
    const volScalarField& TcladOrig, 
    const volScalarField& rhoCoolOrig, 
    const volScalarField& TCoolOrig, 
    const meshToMesh& neutroToFluid)
{

    neutroToFluid.mapTgtToSrc( TfuelOrig, plusEqOp<scalar>(), Tfuel);
    neutroToFluid.mapTgtToSrc( TcladOrig, plusEqOp<scalar>(), Tclad);
    neutroToFluid.mapTgtToSrc( rhoCoolOrig, plusEqOp<scalar>(), rhoCool);
    neutroToFluid.mapTgtToSrc( TCoolOrig, plusEqOp<scalar>(), TCool);

    Tfuel.correctBoundaryConditions();
    Tclad.correctBoundaryConditions();
    rhoCool.correctBoundaryConditions();
    TCool.correctBoundaryConditions();  

}

void Foam::diffusionNeutronics::getFieldsLiquidFuel(
    const volVectorField& UOrig, 
    const volScalarField& porosityOrig, 
    const volScalarField& TfuelOrig, 
    const volScalarField& TcladOrig, 
    const volScalarField& rhoCoolOrig, 
    const volScalarField& TCoolOrig, 
    const meshToMesh& neutroToFluid,
    const rhoThermo& thermo,
    const compressible::turbulenceModel& turb)
{

    neutroToFluid.mapTgtToSrc( TfuelOrig, plusEqOp<scalar>(), Tfuel);
    neutroToFluid.mapTgtToSrc( TcladOrig, plusEqOp<scalar>(), Tclad);
    neutroToFluid.mapTgtToSrc( rhoCoolOrig, plusEqOp<scalar>(), rhoCool);
    neutroToFluid.mapTgtToSrc( TCoolOrig, plusEqOp<scalar>(), TCool);
    neutroToFluid.mapTgtToSrc( UOrig, plusEqOp<vector>(), U);//UNeutro.primitiveFieldRef()
    neutroToFluid.mapTgtToSrc( porosityOrig, plusEqOp<scalar>(), porosity);//porosityNeutro.primitiveFieldRef());

    Tfuel.correctBoundaryConditions();
    Tclad.correctBoundaryConditions();
    rhoCool.correctBoundaryConditions();
    TCool.correctBoundaryConditions();
    U.correctBoundaryConditions();
    porosity.correctBoundaryConditions();

    phi = fvc::flux(U);

    volScalarField diffCoeffOrig = turb.alphat()/thermo.rho()+thermo.mu()/thermo.rho()/xs.ScNo();// (alphaEff=nu/Pr+alphat)
    neutroToFluid.mapTgtToSrc( diffCoeffOrig , plusEqOp<scalar>(), diffCoeffPrec);//.primitiveFieldRef()
    diffCoeffPrec.correctBoundaryConditions();     

}

void Foam::diffusionNeutronics::deformMesh(const meshToMesh& TMToNeutro,const volVectorField& DispOrig)
{

    Info << "Displace neutronic mesh" << endl;

    const volPointInterpolation& neutroMeshPointInterpolation = volPointInterpolation::New(mesh_);

    tmp<pointVectorField> neutroPointsDisplacementOld = neutroMeshPointInterpolation.interpolate(Disp_);

    Disp_*=0.0;
    TMToNeutro.mapSrcToTgt( DispOrig , plusEqOp<vector>(), Disp_.primitiveFieldRef());
    Disp_.correctBoundaryConditions();

    tmp<pointVectorField> neutroPointsDisplacement = neutroMeshPointInterpolation.interpolate(Disp_);

    tmp<pointField> displacedPoints = mesh_.points()
                                    + neutroPointsDisplacement->internalField()
                                    - neutroPointsDisplacementOld->internalField() ;

    mesh_.movePoints(displacedPoints);

}
// ************************************************************************* //

