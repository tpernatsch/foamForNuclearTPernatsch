/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "legacyThermoMechanics.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "coordinateSystem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(legacyThermoMechanics, 0);

    addToRunTimeSelectionTable
    (
        solver,
        legacyThermoMechanics,
        fvMesh
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::legacyThermoMechanics::legacyThermoMechanics
(
    fvMesh& mesh
)
:
    thermoMechanics(mesh),
    solveDisplacement_(this->lookupOrDefault("mechanicsSolver",true)),
    solveThermo_(this->lookupOrDefault("thermalSolver", true)),
    planeStress_
    (
        this->found("rheologyOption") ?
        this->subDict("rheologyOption").get<bool>("planeStress")
        :false
    ),
    linkedFuel_(this->subDict("globalOptions").lookup("linkedFuel")),
    TMEntries_(this->subDict("materials")),
    TMZoneNumber_(TMEntries_.toc().size()),
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
    disp_
    (
        IOobject
        (
            "disp",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    TStruct_
    (
        IOobject
        (
            "TStruct",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TStructRef_
    (
        IOobject
        (
            "TStructRef",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0),
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
        (
            (rhoE_/rho_)/
            (2.0*(1.0 + nu_))
        )*
        twoSymm(fvc::grad(disp_))
    +   (
            nu_*(rhoE_/rho_)/
            (
                (1.0 + nu_)*
                (1.0 - 2.0*nu_)
            )
        )*I*tr(fvc::grad(disp_))
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
    powerDensityNeutronics_
    (
        IOobject
        (
            "powerDensityNeutronics",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    isPorous_ //assume it is not (=0)
    (
        IOobject
        (
            "isPorous",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    nCorr_
    (
        mesh.solutionDict().subDict("stressAnalysis").lookupOrDefault<int>
        (
            "nCorrectors",
            1
        )
    ),
    convergenceTolerance_
    (
        readScalar
        (
            mesh.solutionDict().subDict("stressAnalysis").lookup("D")
        )
    ),
    compactNormalStress_
    (
        mesh.solutionDict().subDict("stressAnalysis").get<bool>
        (
            "compactNormalStress"
        )
    ),
    residual_(0)
{

    PtrList<scalar> TStructRefList(TMZoneNumber_);

    PtrList<scalar> rhoMechList(TMZoneNumber_);
    PtrList<scalar> rhoEList(TMZoneNumber_);
    PtrList<scalar> nuList(TMZoneNumber_);

    PtrList<scalar> CList(TMZoneNumber_);
    PtrList<scalar> rhoKList(TMZoneNumber_);
    PtrList<scalar> alphaList(TMZoneNumber_);

    wordList dictEntries(TMEntries_.toc());


    forAll(dictEntries,zoneI)
    {
        dictionary& dict = TMEntries_.subDict(dictEntries[zoneI]);

        // For compatibility with OFFBEAT material card, just read
        // dimensioned scalars and then assign its value to fields

        dimensionedScalar rho = dict.getOrDefault<dimensionedScalar>("rho",1);
        dimensionedScalar E = dict.getOrDefault<dimensionedScalar>("E", 1);
        dimensionedScalar nu = dict.getOrDefault<dimensionedScalar>("nu", 1);
        dimensionedScalar C = dict.getOrDefault<dimensionedScalar>("Cp", 1);
        dimensionedScalar k = dict.getOrDefault<dimensionedScalar>("k", 1);
        dimensionedScalar alpha = dict.getOrDefault<dimensionedScalar>("alpha", 1);
        dimensionedScalar TStructRef = dict.getOrDefault<dimensionedScalar>("Tref", 0);


        TStructRefList.set(zoneI,new scalar(TStructRef.value()));
        rhoMechList.set(zoneI,new scalar(rho.value()));
        rhoEList.set(zoneI,new scalar(E.value()));
        nuList.set(zoneI,new scalar(nu.value()));
        CList.set(zoneI,new scalar(C.value()));
        rhoKList.set(zoneI,new scalar(k.value()));
        alphaList.set(zoneI,new scalar(alpha.value()));

    }

    // Set volFields based on dictionary
    forAll(dictEntries, zoneI)
    {

        const word& name = dictEntries[zoneI];

        label zoneId = mesh.cellZones().findZoneID(name);

        forAll(mesh.cellZones()[zoneId], cellIlocal)
        {
            label cellIglobal = mesh.cellZones()[zoneId][cellIlocal];

            TStructRef_[cellIglobal] = TStructRefList[zoneI];

            rho_[cellIglobal] = rhoMechList[zoneI];
            rhoE_[cellIglobal] = rhoEList[zoneI];
            nu_[cellIglobal] = nuList[zoneI];

            C_[cellIglobal] = CList[zoneI];
            rhoK_[cellIglobal] = rhoKList[zoneI];
            alpha_[cellIglobal] = alphaList[zoneI];

        }

    }

    TStructRef_.correctBoundaryConditions();
    rho_.correctBoundaryConditions();
    rhoE_.correctBoundaryConditions();
    nu_.correctBoundaryConditions();
    C_.correctBoundaryConditions();
    rhoK_.correctBoundaryConditions();
    alpha_.correctBoundaryConditions();




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

    sigmaD_ =
        mu_*twoSymm(fvc::grad(disp_)) + lambda_*(I*tr(fvc::grad(disp_)));
    divSigmaExp_ = fvc::div(sigmaD_);


    if (compactNormalStress_)
    {
        divSigmaExp_ -=
            fvc::laplacian(2*mu_ + lambda_, disp_, "laplacian(DD,D)");
    }
    else
    {
        divSigmaExp_ -=
            fvc::div((2*mu_ + lambda_)*fvc::grad(disp_), "div(sigmaD)");
    }
    mesh_.setFluxRequired(disp_.name());
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solvers::legacyThermoMechanics::~legacyThermoMechanics()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

// void Foam::solvers::legacyThermoMechanics::getCouplingFieldRefs
// (
//     const objectRegistry& srcTH,
//     const objectRegistry& srcN
// )
// {
//     //- Field names must reflect those defined in createCouplingFields.H
//     TFuelOrig_ =
//         (linkedFuel_) ?
//         srcTH.findObject<volScalarField>("bafflelessTCladAv") :
//         srcTH.findObject<volScalarField>("bafflelessTFuelAv");
//     TStructOrig_ =
//         srcTH.findObject<volScalarField>("bafflelessTStruct");
//     powerDensityOrig_ =
//         srcN.findObject<volScalarField>("powerDensity");
//     //- Initialize mapped fields
//     //this->interpolateCouplingFields(mechToFluid);
// }

// void Foam::solvers::legacyThermoMechanics::interpolateCouplingFields
// (
//     const meshToMesh& mechToFluid,
//     const meshToMesh& mechToNeutro

// )
// {
//     mechToFluid.mapTgtToSrc(*TFuelOrig_, plusEqOp<scalar>(), TFuel_);
//     TFuel_.correctBoundaryConditions();
//     mechToFluid.mapTgtToSrc(*TStructOrig_, plusEqOp<scalar>(), TStructFromTH_);
//     TStructFromTH_.correctBoundaryConditions();

//     mechToNeutro.mapTgtToSrc(*powerDensityOrig_, plusEqOp<scalar>(), powerDensityNeutronics_);
//     powerDensityNeutronics_.correctBoundaryConditions();

// }

void Foam::solvers::legacyThermoMechanics::correctPhysics()
{
    residual_ = 0;
    #include "solveThermalMechanics.H"
}

void Foam::solvers::legacyThermoMechanics::correctTightlyCoupledPhysics()
{
    correctPhysics();
}

scalar Foam::solvers::legacyThermoMechanics::maxDeltaT()
{
    return GREAT;
}




// ************************************************************************* //
