/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2406                                                  |
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
    fvMesh& mesh
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
    initialResidual_(1.0),
    originalPoints_(mesh_.points())
{}


// * * * * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //


autoPtr<fvMesh> Foam::solvers::thermoMechanics::createBaffleLessMesh()
{
   
   
    autoPtr<fvMesh> baffleLessMesh;

    bool removeBaffleBool(false);

    if (mesh_.time().controlDict().found("removeBaffles"))
    {
        const dictionary& removeBafflesDict = mesh_.time().controlDict().subDict("removeBaffles");
       removeBaffleBool = removeBafflesDict.getOrDefault<bool>(mesh_.name(), "false");

    }

    if (removeBaffleBool)
    {
        if (UPstream::nProcs() > 1)
        {
            Info<< "WARNING: the removeBaffles feature is guaranteed to work only "
                << "if each processor contains complete pairs of master and slave "
                << "faces of baffe patches! Otherwise, projection artifacts can "
                << " still occur" << nl << endl;
        }
            
        //- Copy internal geometric data from fluidMesh
        pointField points(mesh_.points());
        faceList faces(mesh_.faces());
        cellList cells(mesh_.cells());

        //- Copy boundary geometric data from fluidMesh
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        wordList patchNames(patches.names());
        List<polyPatch*> pList;
        forAll(patchNames, i)
        {
            pList.append
            (
                new polyPatch
                (
                    patches[patchNames[i]],
                    patches
                )
            );
        }

        //- Assemble mesh as copy from points, faces, cells with no boundary
        baffleLessMesh.reset
        (
            new Foam::fvMesh
            (
                Foam::IOobject
                (
                    IOobject::groupName(mesh_.name(), "baffleLess"),
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::NO_READ
                ),
                std::move(points),
                std::move(faces),
                std::move(cells)
            )
        );

        //- Add boundary data
        baffleLessMesh().addFvPatches(pList);

        //- One of the main issues when passing the actual runTime to removeBaffles
        //  was that the functionObjects execution flags were getting reset to
        //  false. Probably other aspects of runTime were getting modified too. To
        //  avoid all these issues, let removeBaffles operate with a dummy runTime
        //  created on the fly, which is deleted afterwards

        autoPtr<Time> dummyRunTimePtr
        (
            new Time(Foam::Time::controlDictName, mesh_.time().rootPath(), mesh_.time().caseName())
        );
        removeBaffles(baffleLessMesh(), dummyRunTimePtr());
        dummyRunTimePtr.clear();    

        //- Copy mesh zones
        List<pointZone*> pointZonesNB(0);
        List<faceZone*> faceZonesNB(0);
        List<cellZone*> cellZonesNB(0);
        forAll(mesh_.pointZones().names(), i)
        {
            word name(mesh_.pointZones().names()[i]);
            const pointZone& origZone(mesh_.pointZones()[name]);
            pointZonesNB.append
            (
                new pointZone
                (
                    name,
                    origZone,
                    i,
                    mesh_.pointZones()  //- Not sure whether I should use
                                            //  fluidMesh or baffleLessMesh() here,
                                            //  both seem to result in the same
                                            //  zones getting assembled, so, not
                                            //  too big of a deal but this should
                                            //  be checked eventually
                )
            );
        }
        forAll(mesh_.faceZones().names(), i)
        {
            word name(mesh_.faceZones().names()[i]);
            const faceZone& origZone(mesh_.faceZones()[name]);
            faceZonesNB.append
            (
                new faceZone
                (
                    name,
                    origZone,
                    origZone.flipMap(),
                    i,
                    mesh_.faceZones()
                )
            );
        }
        forAll(mesh_.cellZones().names(), i)
        {
            word name(mesh_.cellZones().names()[i]);
            const cellZone& origZone(mesh_.cellZones()[name]);
            cellZonesNB.append
            (
                new cellZone
                (
                    name,
                    origZone,
                    i,
                    mesh_.cellZones()
                )
            );
        }
        baffleLessMesh().addZones(pointZonesNB, faceZonesNB, cellZonesNB);

    }

    else
    {
        baffleLessMesh.set(&mesh_);
    }
    

    return baffleLessMesh;
}

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
                    "multiRegionCouplingDict",
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
                            Info <<regionsFrom[regionFromi]<<" and " << mesh_.name()<<endl;
                            const wordList fieldsList(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("sourceFields")); // list of fields to create
                            const wordList fieldTypes(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("fieldTypes")); // list of types of fields to create
                            fvMesh& baffleLessMesh = const_cast<fvMesh&>(mesh_.time().lookupObject<fvMesh>(mesh_.name()+".baffleLess"));
                            forAll(fieldsList, fieldi)
                            {
                                if (fieldTypes[fieldi] == "scalar")
                                {
                                    volScalarField& field = baffleLessMesh.lookupObjectRef<volScalarField>(fieldsList[fieldi]+".baffleLess");
                                    field.primitiveFieldRef() = mesh_.lookupObject<volScalarField>(fieldsList[fieldi]).primitiveField();
                                    field.correctBoundaryConditions();
                                }
                                else if (fieldTypes[fieldi] == "vector")
                                {
                                    volVectorField& field = baffleLessMesh.lookupObjectRef<volVectorField>(fieldsList[fieldi]+".baffleLess");
                                    field.primitiveFieldRef() = mesh_.lookupObject<volVectorField>(fieldsList[fieldi]).primitiveField();
                                    field.correctBoundaryConditions();
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}



void Foam::solvers::thermoMechanics::deformMesh()
{
    //-Look for the multiRegionDict

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

                tmp<pointVectorField> meshPointsDisplacement = meshPointInterpolation.interpolate(mesh_.lookupObject<volVectorField>(deformDict.subDict(mesh_.name()).get<word>("displacementField")));

                tmp<pointField> displacedPoints =originalPoints_ + meshPointsDisplacement->internalField();

                mesh_.movePoints(displacedPoints);


            }
        }
    }
}



// ************************************************************************* //
