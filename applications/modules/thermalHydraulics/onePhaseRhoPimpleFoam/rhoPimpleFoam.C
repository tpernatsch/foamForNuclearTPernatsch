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

#include "rhoPimpleFoam.H"
#include "addToRunTimeSelectionTable.H"
#include "regimeMapModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(rhoPimpleFoam, 0);
    addToRunTimeSelectionTable
    (
        solver,
        rhoPimpleFoam,
        dynamicFvMesh
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::rhoPimpleFoam::rhoPimpleFoam
(
    dynamicFvMesh& mesh_
)
:
    solver(mesh_),
    mesh_(mesh_),
    pimple_(mesh_),
    pThermo_(fluidThermo::New(mesh_)),
    thermo_(pThermo_()),
    p_
    (
        IOobject
        (
            "p",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    rho_
    (
        IOobject
        (
            "rho",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        thermo_.rho()
    ),
    U_
    (
        IOobject
        (
            "U",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    phi_
    (
        IOobject
        (
            "phi",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho_*U_) & mesh_.Sf()
    ),
    pressureControl_(p_, rho_, pimple_.dict(), false),
    turbulence_
    (
        compressible::turbulenceModel::New
        (
            rho_,
            U_,
            phi_,
            thermo_
        )
    ),
    dpdt_
    (
        IOobject
        (
            "dpdt",
            mesh_.time().timeName(),
            mesh_,
            mesh_.dynamic() ?
                IOobject::READ_IF_PRESENT :
                IOobject::NO_READ,
            mesh_.dynamic() ?
                IOobject::AUTO_WRITE :
                IOobject::NO_WRITE
        ),
        fvc::ddt(p_)
    ),
    K_("K", 0.5*magSqr(U_)),
    MRF_(mesh_),
    rhoMax_("rhoMax", dimDensity, GREAT, pimple_.dict()),
    rhoMin_("rhoMin", dimDensity, Zero, pimple_.dict()),
    fvOptions_(fv::options::New(mesh_)),
    correctPhi_(pimple_.dict().getOrDefault<bool>("correctPhi", false)),
    moveMeshOuterCorrectors_(pimple_.dict().getOrDefault<bool>("moveMeshOuterCorrectors_", false)),
    checkMeshCourantNo_(pimple_.dict().getOrDefault<bool>("checkMeshCourantNo", false)),
    rhoUf_(nullptr),
    residual_(0),
    cumulativeContErr_(0),
    originalPoints_(mesh_.points()),
    solveEnergy_(true),
    stressTensor_
    (
        IOobject
        (
            "stressTensor",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("stressTensor", dimPressure, symmTensor::zero)
    ),
    kappaEff_
    (
        IOobject
        (
            "kappaEff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimPower/dimTemperature/dimLength, 0)
    )
{

    if(mesh_.dynamic())
    {
        Info<< "Constructing face momentum rhoUf_" << endl;
 
        rhoUf_.reset
        (
            new surfaceVectorField
            (
                IOobject
                (
                    "rhoUf_",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fvc::interpolate(rho_*U_)
            )
        );
    }

    // thermo_.validate(args.executable(), "h", "e");
    turbulence_->validate();
    mesh_.setFluxRequired(p_.name());

    solveEnergy_ = pimple_.dict().get<bool>("solveEnergy");

    
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//  Solve according to flags
void Foam::solvers::rhoPimpleFoam::correctPhysics()
{

    residual_=0;
    autoPtr<volScalarField> divrhoU;
    if (correctPhi_)
    {
        divrhoU.reset
        (
            new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi_, rho_, U_))
            )
        );
    }


    // --- Pressure-velocity pimple_ corrector loop
    while (pimple_.loop())
    {
        if (pimple_.firstIter() || moveMeshOuterCorrectors_)
        {
            // Store momentum to set rhoUf_ for introduced faces.
            autoPtr<volVectorField> rhoU;
            if (rhoUf_.valid())
            {
                rhoU.reset(new volVectorField("rhoU", rho_*U_));
            }

            // Do any mesh_ changes -- commented as FSILoop takes care of it
            // mesh_.controlledUpdate();

            if (mesh_.changing())
            {
                MRF_.update();

                if (correctPhi_)
                {
                    // Calculate absolute flux
                    // from the mapped surface velocity
                    phi_ = mesh_.Sf() & rhoUf_();

                    #include "correctPhi.H"

                    // Make the fluxes relative to the mesh_-motion
                    fvc::makeRelative(phi_, rho_, U_);
                }

                // if (checkMeshCourantNo_)
                // {
                //     #include "meshCourantNo.H"
                // }
            }
        }

        if (pimple_.firstIter() && !pimple_.SIMPLErho())
        {
            #include "rhoEqn.H"
        }

        #include "UEqn.H"

        if(solveEnergy_)
            #include "EEqn.H"

        // --- Pressure corrector loop
        while (pimple_.correct())
        {
            if (pimple_.consistent())
            {
                #include "pcEqn.H"
            }
            else
            {
                #include "pEqn.H"
            }
        }

        if (pimple_.turbCorr())
        {
            turbulence_->correct();
        }
    }

    rho_ = thermo_.rho();

    stressTensor_ = -turbulence_->devRhoReff()()-p_*symmTensor(1,0,0,1,0,1);
}

void Foam::solvers::rhoPimpleFoam::correctTightlyCoupledPhysics()
{
    #include "EEqn.H"
}

void Foam::solvers::rhoPimpleFoam::compressibleContinuityErrs()
{
    {
    dimensionedScalar totalMass = fvc::domainIntegrate(rho_);
 
    scalar sumLocalContErr =
        (fvc::domainIntegrate(mag(rho_ - thermo_.rho()))/totalMass).value();
 
    scalar globalContErr =
        (fvc::domainIntegrate(rho_ - thermo_.rho())/totalMass).value();
 
    cumulativeContErr_ += globalContErr;
 
    Info<< "time step continuity errors : sum local = " << sumLocalContErr
        << ", global = " << globalContErr
        << ", cumulative = " << cumulativeContErr_
        << endl;
    }
}


scalar Foam::solvers::rhoPimpleFoam::maxDeltaT()
{

    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
 
    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(phi_))().primitiveField()/rho_.primitiveField()
        );
    
        CoNum = 0.5*gMax(sumPhi/mesh_.V().field())*mesh_.time().deltaTValue();
    
        meanCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*mesh_.time().deltaTValue();
    }
 
    Info<< "Courant Number mean: " << meanCoNum
        << " max: " << CoNum << endl;


    scalar newDeltaT = mesh_.time().controlDict().lookupOrDefault<scalar>("maxDeltaT", GREAT);
    scalar maxCo =mesh_.time().controlDict().lookupOrDefault<scalar>("maxCo", 1.0);
    scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
    scalar deltaTFact =  min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    return min(deltaTFact*mesh_.time().deltaTValue(),newDeltaT);

}

void Foam::solvers::rhoPimpleFoam::correctBaffleLessFields()
{
    if (mesh_.time().controlDict().found("removeBaffles"))
    {
        const dictionary& removeBafflesDict = mesh_.time().controlDict().subDict("removeBaffles");
        if (removeBafflesDict.getOrDefault<bool>(mesh_.name(), false))
        {

            const IOdictionary couplingDict
            (
                IOobject
                (
                    "regionsDict",
                    runTime.time().system(),
                    runTime.db(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                )
            );

            //Lookup for the fields that need to be mapped FROM this mesh
            const dictionary mappingDict(couplingDict.subDict("mappings"));
            //Loop on every region that is not this one and look for the fields in "sourceFields"
            const wordList regions(mappingDict.toc());

            forAll(regions, regioni)
            {
                if (regions[regioni] != mesh_.name()) //look for other regions
                {
                    const dictionary regionFromDict(mappingDict.subDict(regions[regioni]));
                    const wordList regionsFrom(regionFromDict.toc());
                    forAll(regionsFrom, regionFromi)
                    {
                        if (regionsFrom[regionFromi] == mesh_.name())
                        {
                            const wordList fieldsList(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("sourceFields")); // list of fields to create
                            dynamicFvMesh& baffleLessMesh = const_cast<dynamicFvMesh&>(mesh_.time().lookupObject<dynamicFvMesh>(mesh_.name()+".baffleLess"));
                            forAll(fieldsList, fieldi)
                            {
                                correctBaffleLessField<scalar>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<vector>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<tensor>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<symmTensor>(fieldsList[fieldi], baffleLessMesh);
                                correctBaffleLessField<sphericalTensor>(fieldsList[fieldi], baffleLessMesh);

                            }
                        }
                    }
                }
            }
        }
    }
}

template<class Type>
void Foam::solvers::rhoPimpleFoam::correctBaffleLessField(word fieldName, fvMesh& baffleLessMesh)
{

    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if(mesh_.foundObject<VolFieldType>(fieldName))
    {
        VolFieldType& field = baffleLessMesh.lookupObjectRef<VolFieldType>(fieldName+".baffleLess");
        field.primitiveFieldRef() = mesh_.lookupObject<VolFieldType>(fieldName).primitiveField();
        field.correctBoundaryConditions();
    }
}

void Foam::solvers::rhoPimpleFoam::deformMesh()
{
    // Look for the multiRegionDict

    const IOdictionary couplingDict
    (
        IOobject
        (
            "regionsDict",
            runTime.time().system(),
            runTime.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    );

    if (couplingDict.found("meshDeformation"))
    {

        const dictionary deformDict(couplingDict.subDict("meshDeformation"));

        const wordList regions(deformDict.toc());

        forAll(regions, regioni)
        {
            if (regions[regioni] == mesh_.name())
            {
                const volPointInterpolation& meshPointInterpolation = volPointInterpolation::New(mesh_);

                tmp<pointVectorField> meshPointsDisplacement = meshPointInterpolation.interpolate
                (
                    mesh_.lookupObject<volVectorField>
                    (
                        deformDict.subDict(mesh_.name()).get<word>("displacementField")
                    )
                );

                tmp<pointField> displacedPoints = originalPoints_ + meshPointsDisplacement->internalField();

                mesh_.movePoints(displacedPoints);
            }
        }
    }
}



// ************************************************************************* //