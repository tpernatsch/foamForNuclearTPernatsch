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

#include "SNNeutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(SNNeutronics, 0);

    addToRunTimeSelectionTable
    (
        neutronics,
        SNNeutronics,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::SNNeutronics::SNNeutronics
(
    fvMesh& mesh
)
:
    neutronics(mesh),//SNNeutronics is derived from neutronics
    xs_(mesh),
    quadratureSet_
    (
        IOobject
        (
            "quadratureSet",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    discreteDirections_(8*quadratureSet_.lookupOrDefault<int>("discreteDirectionsOct", 1)),
    directionVersors_(discreteDirections_),
    directionWeights_(discreteDirections_),
    legendreMatrices_(xs_.legendreMoments()),
    facePhi_("facePhi",vector(1, 0, 0) & mesh_.Sf()),
    facePhis_(discreteDirections_),
    angularFlux_(xs_.energyGroups()),
    prec_(xs_.precGroups()),
    precStar_(xs_.precGroups()),
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
        dimensionedScalar("", dimensionSet(0,-3,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    flux_(xs_.energyGroups()),
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
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
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
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,-3,-1,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    U_
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
    phi_
    (
        IOobject
        (
            "phi",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvc::flux(U_)
    ),
    porosity_
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
    Tfuel_
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
    Tclad_
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
    rhoCool_
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
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0),0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    diffCoeffPrec_
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
    #include "readQuadratureSet.H"
    #include "calcLegendreMatrices.H"
    #include "createNeutronicsFieldsSN.H"
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::SNNeutronics::~SNNeutronics()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::SNNeutronics::correct(const label couplingIter, scalar& residual, const bool& liquidFuel ) 
{
    #include "solveNeutronicsSN.H"
}


void Foam::SNNeutronics::getFields(
    const volScalarField& TfuelOrig, 
    const volScalarField& TcladOrig, 
    const volScalarField& rhoCoolOrig, 
    const volScalarField& TCoolOrig, 
    const meshToMesh& neutroToFluid)
{

    neutroToFluid.mapTgtToSrc( TfuelOrig, plusEqOp<scalar>(), Tfuel_);
    neutroToFluid.mapTgtToSrc( TcladOrig, plusEqOp<scalar>(), Tclad_);
    neutroToFluid.mapTgtToSrc( rhoCoolOrig, plusEqOp<scalar>(), rhoCool_);
    neutroToFluid.mapTgtToSrc( TCoolOrig, plusEqOp<scalar>(), TCool_);

    Tfuel_.correctBoundaryConditions();
    Tclad_.correctBoundaryConditions();
    rhoCool_.correctBoundaryConditions();
    TCool_.correctBoundaryConditions();  

}

void Foam::SNNeutronics::getFieldsLiquidFuel(
    const volVectorField& UOrig, 
    const volScalarField& porosityOrig, 
    const volScalarField& TfuelOrig, 
    const volScalarField& TcladOrig, 
    const volScalarField& rhoCoolOrig, 
    const volScalarField& TCoolOrig,
    const volScalarField& muOrig, 
    const volScalarField& alphatOrig, 
    const meshToMesh& neutroToFluid
    )
{

    neutroToFluid.mapTgtToSrc( TfuelOrig, plusEqOp<scalar>(), Tfuel_);
    neutroToFluid.mapTgtToSrc( TcladOrig, plusEqOp<scalar>(), Tclad_);
    neutroToFluid.mapTgtToSrc( rhoCoolOrig, plusEqOp<scalar>(), rhoCool_);
    neutroToFluid.mapTgtToSrc( TCoolOrig, plusEqOp<scalar>(), TCool_);
    neutroToFluid.mapTgtToSrc( UOrig, plusEqOp<vector>(), U_);
    neutroToFluid.mapTgtToSrc( porosityOrig, plusEqOp<scalar>(), porosity_);

    Tfuel_.correctBoundaryConditions();
    Tclad_.correctBoundaryConditions();
    rhoCool_.correctBoundaryConditions();
    TCool_.correctBoundaryConditions();
    U_.correctBoundaryConditions();
    porosity_.correctBoundaryConditions();

    phi_ = fvc::flux(U_);

    volScalarField diffCoeffOrig = alphatOrig/rhoCoolOrig+muOrig/rhoCoolOrig/xs_.ScNo();// (alphaEff=nu/Pr+alphat)
    neutroToFluid.mapTgtToSrc( diffCoeffOrig , plusEqOp<scalar>(), diffCoeffPrec_);//.primitiveFieldRef()
    diffCoeffPrec_.correctBoundaryConditions();     

}

// ************************************************************************* //

