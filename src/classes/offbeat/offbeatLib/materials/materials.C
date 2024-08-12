/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

#include "materials.H"
#include "zeroGradientFvPatchField.H"
#include "debug.H"
#include "userParameters.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(materials, 0);
    defineRunTimeSelectionTable(materials, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::materials::materials
(
    const fvMesh& mesh, 
    const dictionary& materialsDict
)
:
    mesh_(mesh),
    materialsDict_(materialsDict),
    T_(nullptr),
    rho_
    (
        IOobject
        (
            "rho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/dimVolume, 0),
        "zeroGradient"
    ),
    Cp_
    (
        IOobject
        (
            "Cp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimSpecificHeatCapacity, 0),
        "zeroGradient"
    ),
    k_
    (
        IOobject
        (
            "k",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPower/dimLength/dimTemperature, 0),
        "zeroGradient"
    ),
    emissivity_
    (
        IOobject
        (
            "emissivity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0),
        "zeroGradient"
    ),
    E_
    (
        IOobject
        (
            "E",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, 0),
        "zeroGradient"
    ),
    nu_
    (
        IOobject
        (
            "nu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, 0),
        "zeroGradient"
    ),
    mu_
    (
        IOobject
        (
            "mu",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, 0),
        "zeroGradient"
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, 0),
        "zeroGradient"
    ),
    threeK_
    (
        IOobject
        (
            "threeK",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, 0),
        "zeroGradient"
    ),
    writeProperties_
    (
        materialsDict_.lookupOrDefault<bool>("writeMaterialProperties",false)
    ),
    maxDeltaMatTol_(GREAT)
    {
        const dictionary& controlDict(mesh_.time().controlDict());


        bool adjustTime = controlDict.lookupOrDefault<bool>("adjustableTimeStep", false);

        if(adjustTime)
        {
            maxDeltaMatTol_ = controlDict.lookupOrDefault<scalar>("maxDeltaMatTol", GREAT);
        }

        if( writeProperties_ )
        {
            // Change write options for all thermomechanical properties fields
            rho_.writeOpt() = IOobject::AUTO_WRITE;
            Cp_.writeOpt() = IOobject::AUTO_WRITE;
            k_.writeOpt() = IOobject::AUTO_WRITE;
            emissivity_.writeOpt() = IOobject::AUTO_WRITE;
            E_.writeOpt() = IOobject::AUTO_WRITE;
            nu_.writeOpt() = IOobject::AUTO_WRITE;
            mu_.writeOpt() = IOobject::AUTO_WRITE;
            lambda_.writeOpt() = IOobject::AUTO_WRITE;
            threeK_.writeOpt() = IOobject::AUTO_WRITE;
        }
    }


Foam::autoPtr<Foam::materials>
Foam::materials::New
(
    const fvMesh& mesh, 
    const dictionary& solverDict
)
{
    // Initialize type for materials class
    word type;

    dictionary materialsDict(solverDict.subOrEmptyDict("materials"));
    solverDict.lookup("materialProperties") >> type;

    Info << "Selecting materials model " << type << endl;

    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("materials::New(const fvMesh&, const dictionary&)")
            << "Unknown materials type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    if (debug)
    {
        Info<< "Selecting materials type "
            << type << endl;
    }

    return autoPtr<materials>(cstrIter()(mesh, materialsDict));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::materials::~materials()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::materials::updateDensity()
{
    if(mesh_.foundObject<fvMesh>("referenceMesh"))
    {
        //- Constant reference to reference volumes
        const fvMesh& refMesh =  mesh_.lookupObject<fvMesh>("referenceMesh");
        const scalarField& Vref =  refMesh.V();

        forAll(mesh_.cellZones(), zoneI)
        {
            const labelList& addr = mesh_.cellZones()[zoneI];

            forAll(addr, i)
            {
                const label cellI = addr[i];

                //- Rescale the cell density with the ref/new volume ratio 
                rho_[cellI] = rho_[cellI] 
                            * Vref[cellI] 
                            / mesh_.V()[cellI];
            }
        }

        rho_.correctBoundaryConditions();
    }
}

Foam::scalar Foam::materials::nextDeltaT()
{
    scalar nextDeltaT(0.0);
    scalar currentDeltaT = mesh_.time().deltaT().value();
    scalar deltaTMultiplier(1.0);

    const scalarField& rhoi = rho_.internalField();
    const scalarField& rhoOld = rho_.oldTime().internalField();
    const scalarField& Cpi = Cp_.internalField();
    const scalarField& CpOld = Cp_.oldTime().internalField();
    const scalarField& ki = k_.internalField();
    const scalarField& k_Old = k_.oldTime().internalField();

    const scalar dRhoPct = max((rhoi - rhoOld) / rhoOld * 100); 
    const scalar dCpPct = max((Cpi - CpOld) / rhoOld * 100);
    const scalar dKPct = max((ki - k_Old) / k_Old * 100);
    
    scalar maxDelta = max(dRhoPct, max(dCpPct, dKPct));
    deltaTMultiplier = maxDeltaMatTol_ / (maxDelta + SMALL);
    nextDeltaT = deltaTMultiplier * currentDeltaT;

    
    Info<< "Maximum deltaT calculated by relative material properties increase: " 
        << mesh_.time().timeToUserTime(nextDeltaT)
        << " with a maximum material properties increase of " << maxDeltaMatTol_*100 << "%" 
        << endl;

    return nextDeltaT;
}

// ************************************************************************* //
