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

#include "neutronics.H"
#include "zeroGradientFvPatchFields.H"
#include "fvmSup.H"
#include "mergeOrSplitBaffles.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{
    defineTypeNameAndDebug(neutronics, 0);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::neutronics::neutronics
(
    fvMesh& mesh
)
:
    solver(mesh),
    IOdictionary
    (
        IOobject
        (
            "neutronicsProperties",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
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
    keff_(reactorState_.lookupOrDefault("keff", 1.0)),
    pTarget_(reactorState_.lookupOrDefault("pTarget", 1.0)),
    pTotOld_(pTarget_),
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
    secondaryPowerDensity_
    (
        IOobject
        (
            "secondaryPowerDensity",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimPower/dimVol, 0.0),
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
        mesh,
        dimensionedScalar("", dimless/dimArea/dimTime, 1),
        zeroGradientFvPatchScalarField::typeName
    ),
    // TFuelOrig_(nullptr),
    // TCladOrig_(nullptr),
    // TCoolOrig_(nullptr),
    // rhoCoolOrig_(nullptr),
    // TStructOrig_(nullptr),
    // TStructMechOrig_(nullptr),
    // UOrig_(nullptr),
    // alphaOrig_(nullptr),
    // alphatOrig_(nullptr),
    // muOrig_(nullptr),
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
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimDensity, SMALL),
        zeroGradientFvPatchScalarField::typeName
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
    TStructMech_
    (
        IOobject
        (
            "TStructMech",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    UPtr_(nullptr),
    alphaPtr_(nullptr),
    alphatPtr_(nullptr),
    muPtr_(nullptr),
    phiPtr_(nullptr),
    diffCoeffPrecPtr_(nullptr),
    // disp_
    // (
    //     IOobject
    //     (
    //         "disp",
    //         mesh.time().timeName(),
    //         mesh,
    //         IOobject::READ_IF_PRESENT,
    //         IOobject::AUTO_WRITE
    //     ),
    //     mesh,
    //     dimensionedVector("d_zero", dimLength, vector::zero),
    //     zeroGradientFvPatchScalarField::typeName
    // ),
    initialResidual_(1.0),
    residual_(0.0),
    eigenvalueNeutronics_
    (
        IOdictionary::lookupOrDefault("eigenvalueNeutronics", false)
    ),
    externalSourceNeutronics_
    (
        IOdictionary::lookupOrDefault("externalSourceNeutronics", false)
    ),
    liquidFuel_
    (
        mesh.time().controlDict().lookupOrDefault("liquidFuel", false)
    ),
    originalPoints_(mesh_.points()),
    externalSource_
    (
        IOobject
        (
            "externalSource",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    externalSourceModulationTimeProfile_
    (
        externalSource_,
        "externalSourceModulationTimeProfile",
        mesh.time()
    )
{
    Info<< "Initial keff = " << keff_ << nl
        << "Is eigenvalue calc : " << eigenvalueNeutronics_ << nl
        << "Is external source calc : " << externalSourceNeutronics_
        << endl;
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::solvers::neutronics::correctBaffleLessFields()
{
    if (mesh_.time().controlDict().found("removeBaffles"))
    {
        const dictionary& removeBafflesDict = mesh_.time().controlDict().subDict("removeBaffles");

        if (removeBafflesDict.get<bool>(mesh_.name()))
        {
            const IOdictionary couplingDict
            (
                IOobject
                (
                    "multiRegionCouplingDict",
                    runTime.time().constant(),
                    runTime.db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            // Lookup for the fields that need to be mapped FROM this mesh
            const dictionary mappingDict(couplingDict.subDict("mappings"));
            // Loop on every region that is not this one and look for the fields in "sourceFields"
            const wordList regions(mappingDict.toc());

            forAll(regions, regioni)
            {
                // Look for other regions
                if (regions[regioni] != mesh_.name())
                {
                    const dictionary regionFromDict(mappingDict.subDict(regions[regioni]));
                    const wordList regionsFrom(regionFromDict.toc());
                    forAll(regionsFrom, regionFromi)
                    {
                        if (regionsFrom[regionFromi] == mesh_.name())
                        {
                            const wordList fieldsList(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("sourceFields")); // list of fields to create
                            fvMesh& baffleLessMesh = const_cast<fvMesh&>(mesh_.time().lookupObject<fvMesh>(mesh_.name()+".baffleLess"));
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
void Foam::solvers::neutronics::correctBaffleLessField(word fieldName, fvMesh& baffleLessMesh)
{

    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if(mesh_.foundObject<VolFieldType>(fieldName))
    {
        VolFieldType& field = baffleLessMesh.lookupObjectRef<VolFieldType>(fieldName+".baffleLess");
        field.primitiveFieldRef() = mesh_.lookupObject<VolFieldType>(fieldName).primitiveField();
        field.correctBoundaryConditions();
    }
}



void Foam::solvers::neutronics::deformMesh()
{
    // Look for the multiRegionDict
    const IOdictionary couplingDict
    (
        IOobject
        (
            "multiRegionCouplingDict",
            runTime.time().constant(),
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

                tmp<pointVectorField> meshPointsDisplacement
                (
                    meshPointInterpolation.interpolate
                    (
                        mesh_.lookupObject<volVectorField>
                        (
                            deformDict.subDict(mesh_.name()).get<word>("displacementField")
                        )
                    )
                );

                tmp<pointField> displacedPoints = originalPoints_ + meshPointsDisplacement->internalField();

                mesh_.movePoints(displacedPoints);

                Info<< "Deforming " << mesh_.name()
                    << " mesh according to displacement field "
                    << deformDict.subDict(mesh_.name()).get<word>("displacementField")
                    << endl;
            }
        }
    }
}


// ************************************************************************* //
