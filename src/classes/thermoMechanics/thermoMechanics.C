/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
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

#include "thermoMechanics.H"
#include "zeroGradientFvPatchFields.H"
#include "mergeOrSplitBaffles.H"
#include "fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace solvers
{

    defineTypeNameAndDebug(thermoMechanics, 0);
    // defineRunTimeSelectionTable(thermoMechanics, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solvers::thermoMechanics::thermoMechanics
(
    dynamicFvMesh& mesh
)
:
    solver(mesh),
    IOdictionary
    (
        IOobject
        (
            "thermoMechanicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    meshDisp_
    (
        IOobject
        (
            "meshDisp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("", dimLength, vector::zero),
        zeroGradientFvPatchScalarField::typeName
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
    fuelOrientation_(this->subDict("globalOptions").lookup("pinDirection")),
    TStructFromTH_
    (
        IOobject
        (
            "TStructFromTH",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TFuel_
    (
        IOobject
        (
            "TFuel",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("", dimensionSet(0,0,0,1,0,0,0), 1.0),
        zeroGradientFvPatchScalarField::typeName
    ),
    TFuelRef_
    (
        IOobject
        (
            "TFuelRef",
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
    TCRRef_
    (
        IOobject
        (
            "TCRRef",
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
    TMEntries_(this->subDict("materials")),
    zonesNumber_(TMEntries_.toc().size()),
    initialResidual_(1.0),
    originalPoints_(mesh_.points())

{
    PtrList<scalar> TFuelRefList(zonesNumber_);
    PtrList<scalar> alphaFuelList(zonesNumber_);
    PtrList<scalar> TCRRefList(zonesNumber_);
    PtrList<scalar> alphaCRList(zonesNumber_);

    wordList dictEntries(TMEntries_.toc());


    forAll(dictEntries,zoneI)
    {
        dictionary& dict = TMEntries_.subDict(dictEntries[zoneI]);

        dimensionedScalar TFuelRef = dict.getOrDefault<dimensionedScalar>("TFuelRef", 0.0);
        dimensionedScalar alphaFuel = dict.getOrDefault<dimensionedScalar>("alphaFuel", 0.0);
        dimensionedScalar TCRRef = dict.getOrDefault<dimensionedScalar>("TCRRef", 0.0);
        dimensionedScalar alphaCR = dict.getOrDefault<dimensionedScalar>("alphaCR", 0.0);


        TFuelRefList.set
        (
            zoneI,
            new scalar(TFuelRef.value())
        );
        alphaFuelList.set
        (
            zoneI,
            new scalar(alphaFuel.value())
        );

        TCRRefList.set(zoneI,new scalar(TCRRef.value()));
        alphaCRList.set(zoneI,new scalar(alphaCR.value()));
    }

    // Set volFields based on dictionary
    forAll(dictEntries, zoneI)
    {

        const word& name = dictEntries[zoneI];

        label zoneId = mesh.cellZones().findZoneID(name);

        forAll(mesh.cellZones()[zoneId], cellIlocal)
        {
            label cellIglobal = mesh.cellZones()[zoneId][cellIlocal];

            TFuelRef_[cellIglobal] = TFuelRefList[zoneI];
            alphaFuel_[cellIglobal] = alphaFuelList[zoneI];
            TCRRef_[cellIglobal] = TCRRefList[zoneI];
            alphaCR_[cellIglobal] = alphaCRList[zoneI];
        }

    }

    TFuelRef_.correctBoundaryConditions();
    alphaFuel_.correctBoundaryConditions();

    TCRRef_.correctBoundaryConditions();
    alphaCR_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

void Foam::solvers::thermoMechanics::correctBaffleLessFields()
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
                    "regionsDict",
                    runTime.time().constant(),
                    runTime.db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );

            //Lookup for the fields that need to be mapped FROM this mesh
            const dictionary mappingDict(couplingDict.subDict("mappings"));
            //Loop on every region that is not this one and look for the fields in "sourceFields"
            const wordList regions(mappingDict.toc());

            forAll(regions, regioni)
            {
                if(regions[regioni]!=mesh_.name()) //look for other regions
                {
                    const dictionary regionFromDict(mappingDict.subDict(regions[regioni]));
                    const wordList regionsFrom(regionFromDict.toc());
                    forAll(regionsFrom, regionFromi)
                    {
                        if(regionsFrom[regionFromi]==mesh_.name())
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
void Foam::solvers::thermoMechanics::correctBaffleLessField(word fieldName, fvMesh& baffleLessMesh)
{

    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    if(mesh_.foundObject<VolFieldType>(fieldName))
    {
        VolFieldType& field = baffleLessMesh.lookupObjectRef<VolFieldType>(fieldName+".baffleLess");
        field.primitiveFieldRef() = mesh_.lookupObject<VolFieldType>(fieldName).primitiveField();
        field.correctBoundaryConditions();
    }
}



void Foam::solvers::thermoMechanics::deformMesh()
{
    //-Look for the multiRegionDict

    const IOdictionary couplingDict
    (
        IOobject
        (
            "regionsDict",
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

                tmp<pointVectorField> meshPointsDisplacement = meshPointInterpolation.interpolate(mesh_.lookupObject<volVectorField>(deformDict.subDict(mesh_.name()).get<word>("displacementField")));

                tmp<pointField> displacedPoints =originalPoints_ + meshPointsDisplacement->internalField();

                mesh_.movePoints(displacedPoints);


            }
        }
    }
}



// ************************************************************************* //
