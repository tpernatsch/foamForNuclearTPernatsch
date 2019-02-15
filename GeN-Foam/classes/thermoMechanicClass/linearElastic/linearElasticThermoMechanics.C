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

#include "linearElasticThermoMechanics.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearElasticThermoMechanics, 0);

    addToRunTimeSelectionTable
    (
        thermoMechanics,
        linearElasticThermoMechanics,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::linearElasticThermoMechanics::linearElasticThermoMechanics
(
    fvMesh& mesh
)
:
    thermoMechanics(mesh),
    planeStress_(this->lookup("planeStress")),
    linkedFuel_(this->lookup("linkedFuel")),
    fuelOrientation_(this->lookup("fuelOrientation")),
    TrefStructures_("",dimensionSet(0,0,0,1,0,0,0),this->lookup("TrefStructure")),
    TMEntries_(this->lookup("zones")),
    TMZoneNumber_(TMEntries_.size()),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(1,-3,0,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoE_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(1,-1,-2,0,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    nu_
    (
        IOobject
        (
            "nu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    threeKalpha_
    (
        IOobject
        (
            "threeKalpha",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    DT_
    (
        IOobject
        (
            "DT",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 2, -1 , 0, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    C_
    (
        IOobject
        (
            "C",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 2, -2 , -1, 0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    rhoK_
    (
        IOobject
        (
            "rhoK",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(1, 1, -3 , -1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alpha_
    (
        IOobject
        (
            "alpha",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0 , -1, 0), 1.0e-5),
        zeroGradientFvPatchScalarField::typeName
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,2,-2,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),   
    mu_
    (
        IOobject
        (
            "mu",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,2,-2,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,2,-2,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    threeK_
    (
        IOobject
        (
            "threeK",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,2,-2,0,0,0,0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TrefFuel_
    (
        IOobject
        (
            "TrefFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0 , 1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaFuel_
    (
        IOobject
        (
            "alphaFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0 , -1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TrefCR_
    (
        IOobject
        (
            "TrefCR",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0 , 1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    alphaCR_
    (
        IOobject
        (
            "alphaCR",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("0", dimensionSet(0, 0, 0 , -1, 0), 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    Disp_
    (
        IOobject
        (
            "Disp",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    fuelDisp_
    (
        IOobject
        (
            "fuelDisp",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    CRDisp_
    (
        IOobject
        (
            "CRDisp",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    fuelDispVector_
    (
        IOobject
        (
            "fuelDispVector",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        (fuelDisp_*vector(0,0,1))
    ),
    T_
    (
        IOobject
        (
            "T",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 1.0),//TrefStructures,//important to impose correct BC (this value will remain in the soft region surrounding everything)
        zeroGradientFvPatchScalarField::typeName
    ),
    sigmaD_
    (
        IOobject
        (
            "sigmaD",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        ((rhoE_/rho_)/(2.0*(1.0 + nu_)))*twoSymm(fvc::grad(Disp_)) + (nu_*(rhoE_/rho_)/((1.0 + nu_)*(1.0 - 2.0*nu_)))*(I*tr(fvc::grad(Disp_)))
    ),
    divSigmaExp_
    (
        IOobject
        (
            "divSigmaExp",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::div(sigmaD_)
    ),
    TavFuel_
    (
        IOobject
        (
            "TavFuelMech",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    gapWidth
    (
        IOobject
        (
            "gapWidth",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("gapWidth", dimensionSet(0,1,0,0,0,0,0), 0.0),
        calculatedFvPatchField<scalar>::typeName
    ),
    nCorr_(mesh.solutionDict().subDict("stressAnalysis").lookupOrDefault<int>("nCorrectors", 1)),
    convergenceTolerance_(readScalar(mesh.solutionDict().subDict("stressAnalysis").lookup("D"))),
    compactNormalStress_(mesh.solutionDict().subDict("stressAnalysis").lookup("compactNormalStress"))
{

    PtrList<scalar > rhoMechList(TMZoneNumber_);
    PtrList<scalar > rhoEList(TMZoneNumber_);
    PtrList<scalar > nuList(TMZoneNumber_);

    PtrList<scalar > CList(TMZoneNumber_);
    PtrList<scalar > rhoKList(TMZoneNumber_);
    PtrList<scalar > alphaList(TMZoneNumber_);

    PtrList<scalar > TrefFuelList(TMZoneNumber_);
    PtrList<scalar > alphaFuelList(TMZoneNumber_);

    PtrList<scalar > TrefCRList(TMZoneNumber_);
    PtrList<scalar > alphaCRList(TMZoneNumber_);

    forAll(TMEntries_,zoneI)
    {
        dictionary& dict = TMEntries_[zoneI].dict();

        rhoMechList.set(zoneI,new scalar(dict.lookupOrDefault("rho",1.0)));
        rhoEList.set(zoneI,new scalar(dict.lookupOrDefault("E",1.0)));
        nuList.set(zoneI,new scalar(dict.lookupOrDefault("nu",0.0)));

        CList.set(zoneI,new scalar(dict.lookupOrDefault("C",1.0)));
        rhoKList.set(zoneI,new scalar(dict.lookupOrDefault("k",1.0)));
        alphaList.set(zoneI,new scalar(dict.lookupOrDefault("alpha",0.0)));

        TrefFuelList.set(zoneI,new scalar(dict.lookupOrDefault("TrefFuel",0.0)));
        alphaFuelList.set(zoneI,new scalar(dict.lookupOrDefault("alphaFuel",0.0)));

        TrefCRList.set(zoneI,new scalar(dict.lookupOrDefault("TrefCR",0.0)));
        alphaCRList.set(zoneI,new scalar(dict.lookupOrDefault("alphaCR",0.0)));
    }


    // Set volFields based on dictionary
    forAll(TMEntries_, zoneI)
    {

        const word& name = TMEntries_[zoneI].keyword();

        label zoneId = mesh.cellZones().findZoneID(name);

        forAll(mesh.cellZones()[zoneId], cellIlocal)
        {
            label cellIglobal = mesh.cellZones()[zoneId][cellIlocal];

            rho_[cellIglobal] = rhoMechList[zoneI];
            rhoE_[cellIglobal] = rhoEList[zoneI];
            nu_[cellIglobal] = nuList[zoneI];

            C_[cellIglobal] = CList[zoneI];
            rhoK_[cellIglobal] = rhoKList[zoneI];
            alpha_[cellIglobal] = alphaList[zoneI];

            TrefFuel_[cellIglobal] = TrefFuelList[zoneI];
            alphaFuel_[cellIglobal] = alphaFuelList[zoneI];

            TrefCR_[cellIglobal] = TrefCRList[zoneI];
            alphaCR_[cellIglobal] = alphaCRList[zoneI];
        }

    }

    rho_.correctBoundaryConditions();
    rhoE_.correctBoundaryConditions();
    nu_.correctBoundaryConditions();

    C_.correctBoundaryConditions();
    rhoK_.correctBoundaryConditions();
    alpha_.correctBoundaryConditions();

    TrefFuel_.correctBoundaryConditions();
    alphaFuel_.correctBoundaryConditions();

    TrefCR_.correctBoundaryConditions();
    alphaCR_.correctBoundaryConditions();


    E_ = rhoE_/rho_ ;
    mu_ = E_/(2.0*(1.0 + nu_)) ;
    lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - 2.0*nu_)) ;
    threeK_ = E_/(1.0 - 2.0*nu_) ;


    mu_.correctBoundaryConditions();
    lambda_.correctBoundaryConditions();
    threeK_.correctBoundaryConditions();
    sigmaD_.correctBoundaryConditions();

    if (planeStress_)
    {
        Info<< "Plane Stress\n" << endl;

        lambda_ = nu_*E_/((1.0 + nu_)*(1.0 - nu_));
        threeK_ = E_/(1.0 - nu_);
        lambda_.correctBoundaryConditions();
        threeK_.correctBoundaryConditions();        
    }
    else
    {
        Info<< "Plane Strain\n" << endl;
    }
   

    Info<< "Normalising k : k/rho\n" << endl;
    volScalarField k(rhoK_/rho_);

    Info<< "Calculating thermal coefficients\n" << endl;

    threeKalpha_ = threeK_*alpha_;
    DT_ = k/C_;

    sigmaD_ =   mu_*twoSymm(fvc::grad(Disp_)) + lambda_*(I*tr(fvc::grad(Disp_)));
    divSigmaExp_ = fvc::div(sigmaD_);


    if (compactNormalStress_)
    {
        divSigmaExp_ -= fvc::laplacian(2*mu_ + lambda_, Disp_, "laplacian(DD,D)");
    }
    else
    {
        divSigmaExp_ -= fvc::div((2*mu_ + lambda_)*fvc::grad(Disp_), "div(sigmaD)");
    }
    mesh_.setFluxRequired(Disp_.name());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearElasticThermoMechanics::~linearElasticThermoMechanics()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::linearElasticThermoMechanics::correct(scalar& residual) 
{
    #include "solveThermalMechanics.H"
}

void Foam::linearElasticThermoMechanics::getFields
(
const volScalarField& T, 
const volScalarField& TavFuel,
const volScalarField& TavClad,
const meshToMesh& TMToFluid) 
{

    TMToFluid.mapTgtToSrc( T, plusEqOp<scalar>(), T_);


    TMToFluid.mapTgtToSrc(T, plusEqOp<scalar>(), T_);
    if(linkedFuel_)
    {
        TMToFluid.mapTgtToSrc(TavClad, plusEqOp<scalar>(), TavFuel_);
    }
    else
    {
        TMToFluid.mapTgtToSrc(TavFuel, plusEqOp<scalar>(), TavFuel_);
    }

    T_.correctBoundaryConditions();
    TavFuel_.correctBoundaryConditions();

    T_.correctBoundaryConditions();
    TavFuel_.correctBoundaryConditions();
}


// ************************************************************************* //

