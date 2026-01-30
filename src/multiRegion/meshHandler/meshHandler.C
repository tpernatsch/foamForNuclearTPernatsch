/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
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

#include "meshHandler.H"
#include "solver.H"
#include "Time.H"
#include "interpolationCellPoint.H"
#include "radialBasisFunctionInterpolation.H"
#include "mergeOrSplitBaffles.H"
#include "hexCellFvMesh.H"
#include "dynamicFvMesh.H"
#include "staticFvMesh.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshHandler::meshHandler(const Time& runTime)
{
    // Read meshHandler dict. Will be used later
    IOdictionary multiRegionCouplingDict
    (
        IOobject
        (
            "regionsDict",
            runTime.time().system(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    IOdictionary controlDict
    (
        IOobject
        (
            "controlDict",
            runTime.system(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    mappingDict_ = multiRegionCouplingDict.subDict("mappings");

    // First create all the meshes -> they are created and handled by meshHandler

    wordList allRegions(multiRegionCouplingDict.subDict("regionSolvers").subDict("Level_0").toc());
    wordList meshNames(0);

    wordList availableSolvers =
    (
        Foam::solver::dynamicFvMeshConstructorTablePtr_->toc()
    );   

        static const HashSet<word> loopTypes = { "picardLoop", "picardLoopNoFluid", "FSILoop", "CHTLoop", "multiScaleLoop"};

    forAll(allRegions, regioni)
    {
        word solverName(multiRegionCouplingDict.subDict("regionSolvers").subDict("Level_0").get<word>(allRegions[regioni]));
        if (!(loopTypes.found(solverName)))
        {
            if(availableSolvers.found(solverName))
            {
                meshNames.append(allRegions[regioni]);
            }
            else
            {
                FatalErrorInFunction <<
                "Solver " << solverName <<
                " does not exist. Possible solvers are" <<
                Foam::solver::dynamicFvMeshConstructorTablePtr_->toc()
                <<endl<<exit(FatalError);
            }
        }
    }

    // If multiPhysicsSolvers are present, also add those extra solvers

    forAll(multiRegionCouplingDict.subDict("regionSolvers").toc(), MPsolvI) // loop over all multiphysics solvers
    {
        if(multiRegionCouplingDict.subDict("regionSolvers").toc()[MPsolvI] != "Level_0")
        {
            wordList subSolvers // get list of solvers belonging to multiphysics solver i
            (
                multiRegionCouplingDict
                    .subDict("regionSolvers")
                    .subDict(multiRegionCouplingDict.subDict("regionSolvers").toc()[MPsolvI])
                    .subDict("subSolvers").toc()
            );
    
            forAll(subSolvers, solvI) // loop over subsolvers and add them to list of meshes to create only if they are not sub-multiphysicssolvers
            {
                if (!(loopTypes.found(multiRegionCouplingDict
                    .subDict("regionSolvers")
                    .subDict(multiRegionCouplingDict.subDict("regionSolvers").toc()[MPsolvI])
                    .subDict("subSolvers")
                    .get<word>(subSolvers[solvI])))
                )
                {
                    meshNames.append(subSolvers[solvI]);


                    // Check if the solver is a multiScaleLoop

                    if(multiRegionCouplingDict
                       .subDict("regionSolvers")
                       .subDict(multiRegionCouplingDict.subDict("regionSolvers").toc()[MPsolvI])
                       .found("multiScaleLoopCoeffs")
                    )
                    {
                        dictionary mSLCDict = 
                        (
                            multiRegionCouplingDict
                            .subDict("regionSolvers")
                            .subDict(multiRegionCouplingDict.subDict("regionSolvers").toc()[MPsolvI])
                            .subDict("multiScaleLoopCoeffs")
                        );

                        vectorList locations = mSLCDict.get<vectorList>("locations");

                        word regionToReplicate = mSLCDict.get<word>("regionToReplicate");

                        forAll(locations, locI)
                        {
                            meshNames.append(regionToReplicate + Foam::name(locI));
                        }

                    }

                }
                else
                {
                    word solverIName(multiRegionCouplingDict
                    .subDict("regionSolvers")
                    .subDict(multiRegionCouplingDict.subDict("regionSolvers").toc()[MPsolvI])
                    .subDict("subSolvers")
                    .get<word>(subSolvers[solvI]));

                    if(!(availableSolvers.found(solverIName)))
                    {
                        FatalErrorInFunction <<
                        "Solver " << solverIName <<
                        " does not exist. Possible solvers are" <<
                        Foam::solver::dynamicFvMeshConstructorTablePtr_->toc()
                        <<endl<<exit(FatalError);
                    }
                }
            }
        }
    }

    meshes_.setSize(meshNames.size());

    forAll(meshNames, meshI)
    {
        meshes_.set
        (
            meshI,
            Foam::dynamicFvMesh::New
            (
                IOobject
                (
                    meshNames[meshI],
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }

    // Set size of list of (possibly) baffleLess meshes
    mappingMeshes_.setSize(meshes_.size());
    // Set size of list of couplingFields
    scalarCouplingFields_.setSize(meshes_.size());
    vectorCouplingFields_.setSize(meshes_.size());
    couplingFields_.setSize(meshes_.size());

    createBaffleLessMeshes(runTime);
    //createCouplingFields(runTime);


    // Create mappings
    mappingList_.setSize(meshes_.size());

    forAll(meshes_, i)
    {
        mappingList_.set
        (
            i,
            new PtrList<meshToMesh>
            (
                meshes_.size()
            )
        );
    }


    // Loop over all regions
    forAll(meshes_, i)
    {
        // Check if i find entry in the meshHandler dict
        if (mappingDict_.found(meshes_[i].name()))
        {
            // Loop over all other regions
            forAll(meshes_, j)
            {
                // Check if i find entry in the meshHandler subDict
                if (mappingDict_.subDict(meshes_[i].name()).found(meshes_[j].name()))
                {
                    bool removeBaffles(false);
                    if(runTime.controlDict().found("removeBaffles"))
                    {
                        removeBaffles= runTime.controlDict().
                        subDict("removeBaffles").
                        getOrDefault<bool>(meshes_[j].name(),false);
                    }
                    mappingList_[i].set
                    (
                        j,
                        new meshToMesh
                        (
                            meshes_[i],
                            removeBaffles
                                ? mappingMeshes_[j]
                                : meshes_[j],
                            Foam::meshToMesh::interpolationMethod::imCellVolumeWeight, // for now hard-coded, possibly implement meshHandler type via dict
                            Foam::meshToMesh::procMapMethod::pmAABB,
                            false
                        )
                    );
                }
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshHandler::~meshHandler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::meshHandler::mapAllFields(const Time& runTime)
{
    forAll(meshes_, i)
    {
        forAll(meshes_, j)
        {
            if (mappingDict_.subDict(meshes_[i].name()).found(meshes_[j].name()))
            {
                dictionary iTojRegionmappingDict_(mappingDict_.subDict(meshes_[i].name()).subDict(meshes_[j].name()));
                List<word> targetFields(iTojRegionmappingDict_.getOrDefault<List<word>>("targetFields", List<word>::null()));
                List<word> sourceFields(iTojRegionmappingDict_.getOrDefault<List<word>>("sourceFields",List<word>::null()));
                // List<word> fieldTypes(iTojRegionmappingDict_.getOrDefault<List<word>>("fieldTypes",List<word>::null()));

                bool removeBaffles(false);

                if (runTime.controlDict().found("removeBaffles"))
                {
                    const dictionary& removeBafflesDict(runTime.controlDict().subDict("removeBaffles"));

                    if (removeBafflesDict.getOrDefault<bool>(meshes_[j].name(), false))
                    {
                        forAll(sourceFields, sourcei)
                        {
                            sourceFields[sourcei]+=".baffleLess";
                        }
                        removeBaffles = true;
                    }
                }

                if (targetFields.size()!=sourceFields.size())
                {
                    FatalErrorInFunction
                        << "Number of target fields from region "
                        << meshes_[j].name()
                        << " is not equal to number of source fields to region "
                        << mappingMeshes_[i].name() << " !"
                        << exit(FatalError);
                }
                else
                {
                    forAll(sourceFields, fieldi)
                    {
                        map<scalar>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                        map<vector>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                        map<tensor>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                        map<symmTensor>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                    }
                }
            }
        }
    }

}


void Foam::meshHandler::mapTheseFields(const Time& runTime, wordList meshToMap)
{
    forAll(meshes_, i)
    {
        forAll(meshes_, j)
        {
            if(meshToMap.contains(meshes_[i].name()) and meshToMap.contains(meshes_[j].name()))
            {
                if(mappingDict_.subDict(meshes_[i].name()).found(meshes_[j].name()))
                {
                    dictionary iTojRegionmappingDict_(mappingDict_.subDict(meshes_[i].name()).subDict(meshes_[j].name()));
                    List<word> targetFields(iTojRegionmappingDict_.getOrDefault<List<word>>("targetFields", List<word>::null()));
                    List<word> sourceFields(iTojRegionmappingDict_.getOrDefault<List<word>>("sourceFields",List<word>::null()));

                    bool removeBaffles(false);

                    if (runTime.controlDict().found("removeBaffles"))
                    {
                        const dictionary& removeBafflesDict(runTime.controlDict().subDict("removeBaffles"));

                        if (removeBafflesDict.getOrDefault<bool>(meshes_[j].name(), false))
                        {
                            forAll(sourceFields, sourcei)
                            {
                                sourceFields[sourcei]+=".baffleLess";
                            }
                            removeBaffles = true;
                        }
                    }

                    if (targetFields.size() != sourceFields.size())
                    {
                        FatalErrorInFunction
                            << "Number of target fields from region "
                            << meshes_[j].name()
                            << " is not equal to number of source fields to region "
                            << mappingMeshes_[i].name() << "!"
                            << exit(FatalError);
                    }
                    else
                    {
                        forAll(sourceFields, fieldi)
                        {
                            map<scalar>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                            map<vector>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                            map<tensor>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                            map<symmTensor>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles);
                        }
                    }
                }
            }
        }
    }

}



void Foam::meshHandler::initializeMappedFields
(
    const Time& runTime,
    const bool initializeMappedFields
)
{
    forAll(meshes_, i)
    {
        forAll(meshes_, j)
        {
            if( mappingDict_.subDict(meshes_[i].name()).found(meshes_[j].name()))
            {
                dictionary iTojRegionmappingDict_(mappingDict_.subDict(meshes_[i].name()).subDict(meshes_[j].name()));
                List<word> targetFields(iTojRegionmappingDict_.getOrDefault<List<word>>("targetFields", List<word>::null()));
                List<word> sourceFields(iTojRegionmappingDict_.getOrDefault<List<word>>("sourceFields",List<word>::null()));

                bool removeBaffles(false);

                if (runTime.controlDict().found("removeBaffles"))
                {
                    const dictionary& removeBafflesDict(runTime.controlDict().subDict("removeBaffles"));

                    if (removeBafflesDict.getOrDefault<bool>(meshes_[j].name(), false))
                    {
                        forAll(sourceFields, sourcei)
                        {
                            sourceFields[sourcei] += ".baffleLess";
                        }
                        removeBaffles = true;
                    }
                }

                if (targetFields.size() != sourceFields.size())
                {
                    FatalErrorInFunction
                        << "Number of target fields from region "
                        << meshes_[j].name()
                        << " is not equal to number of source fields to region "
                        << mappingMeshes_[i].name()<<"!"
                        << exit(FatalError);
                }
                else
                {
                    forAll(sourceFields, fieldi)
                    {
                        mapAndWrite<scalar>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles, initializeMappedFields);
                        mapAndWrite<vector>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles, initializeMappedFields);
                        mapAndWrite<tensor>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles, initializeMappedFields);
                        mapAndWrite<symmTensor>(sourceFields[fieldi], targetFields[fieldi],i, j, removeBaffles, initializeMappedFields);
                    }
                }
            }
        }
    }
}


void Foam::meshHandler::createBaffleLessMeshes(const Time& runTime)
{
    forAll(meshes_, regioni)
    {
        autoPtr<fvMesh> baffleLessMesh;

        bool removeBaffleBool(false);

        if (runTime.controlDict().found("removeBaffles"))
        {
            const dictionary& removeBafflesDict = runTime.controlDict().subDict("removeBaffles");
            removeBaffleBool = removeBafflesDict.getOrDefault<bool>(meshes_[regioni].name(), "false");
        }

        if (removeBaffleBool)
        {

            Info<< "Creating baffleless mesh for region " << meshes_[regioni].name()
                << endl;
            if (UPstream::nProcs() > 1)
            {
                Info<< "WARNING: the removeBaffles feature is guaranteed to work only "
                    << "if each processor contains complete pairs of master and slave "
                    << "faces of baffe patches! Otherwise, projection artifacts can "
                    << " still occur" << nl << endl;
            }


            //  Copy internal geometric data from fluidMesh
            pointField points(meshes_[regioni].points());
            faceList faces(meshes_[regioni].faces());
            cellList cells(meshes_[regioni].cells());

            //  Copy boundary geometric data from fluidMesh
            const polyBoundaryMesh& patches = meshes_[regioni].boundaryMesh();
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

            //  Assemble mesh as copy from points, faces, cells with no boundary
            baffleLessMesh.reset
            (
                new Foam::fvMesh
                (
                    Foam::IOobject
                    (
                        IOobject::groupName(meshes_[regioni].name(), "baffleLess"),
                        runTime.timeName(),
                        runTime,
                        Foam::IOobject::NO_READ
                    ),
                    std::move(points),
                    std::move(faces),
                    std::move(cells)
                )
            );

            

            //  Add boundary data
            baffleLessMesh().addFvPatches(pList);

            //  One of the main issues when passing the actual runTime to removeBaffles
            //  was that the functionObjects execution flags were getting reset to
            //  false. Probably other aspects of runTime were getting modified too. To
            //  avoid all these issues, let removeBaffles operate with a dummy runTime
            //  created on the fly, which is deleted afterwards

            autoPtr<Time> dummyRunTimePtr
            (
                new Time
                (
                    Foam::Time::controlDictName,
                    meshes_[regioni].time().rootPath(),
                    meshes_[regioni].time().caseName()
                )
            );
            removeBaffles(baffleLessMesh(), dummyRunTimePtr());
            dummyRunTimePtr.clear();

            //  Copy mesh zones
            List<pointZone*> pointZonesNB(0);
            List<faceZone*> faceZonesNB(0);
            List<cellZone*> cellZonesNB(0);
            forAll(meshes_[regioni].pointZones().names(), i)
            {
                word name(meshes_[regioni].pointZones().names()[i]);
                const pointZone& origZone(meshes_[regioni].pointZones()[name]);
                pointZonesNB.append
                (
                    new pointZone
                    (
                        name,
                        origZone,
                        i,
                        meshes_[regioni].pointZones()
                    )
                );
            }
            forAll(meshes_[regioni].faceZones().names(), i)
            {
                word name(meshes_[regioni].faceZones().names()[i]);
                const faceZone& origZone(meshes_[regioni].faceZones()[name]);
                faceZonesNB.append
                (
                    new faceZone
                    (
                        name,
                        origZone,
                        origZone.flipMap(),
                        i,
                        meshes_[regioni].faceZones()
                    )
                );
            }
            forAll(meshes_[regioni].cellZones().names(), i)
            {
                word name(meshes_[regioni].cellZones().names()[i]);
                const cellZone& origZone(meshes_[regioni].cellZones()[name]);
                cellZonesNB.append
                (
                    new cellZone
                    (
                        name,
                        origZone,
                        i,
                        meshes_[regioni].cellZones()
                    )
                );
            }
            baffleLessMesh().addZones(pointZonesNB, faceZonesNB, cellZonesNB);
        }
        else
        {
            baffleLessMesh = nullptr;
        }

        mappingMeshes_.set(regioni,baffleLessMesh);
    }
}


void Foam::meshHandler::createCouplingFields(const Time& runTime)
{

    const IOdictionary couplingDict
    (
        IOobject
        (
            "regionsDict",
            runTime.time().system(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const dictionary mappingDict(couplingDict.subDict("mappings"));
    const wordList regionsInDict(mappingDict.toc());

    forAll(meshes_, regioni) //iterate over the meshes
    {
        if (runTime.controlDict().found("removeBaffles"))
        {
            const dictionary& removeBafflesDict = runTime.controlDict().subDict("removeBaffles");
            if (removeBafflesDict.getOrDefault<bool>(meshes_[regioni].name(), false))
            {

                label nCouplingFields(0);

                couplingFields_[regioni].setSize(nCouplingFields+1);

                forAll(regionsInDict, i) //iterate over the regions found in the dictionary
                {
                    if (regionsInDict[i]!=meshes_[regioni].name()) // look for other (than regionMehses_[regioni]) regions source fields
                    {
                        const dictionary regionFromDict(mappingDict.subDict(regionsInDict[i]));
                        const wordList regionsFrom(regionFromDict.toc());
                        // scalarCouplingFields_[regioni].setSize(regionsFrom.size());

                        forAll(regionsFrom, regionFromi)// loop over regions in the sub-dict
                        {
                            if(regionsFrom[regionFromi] == meshes_[regioni].name())// if i find regionMehses_[regioni] in the subdict
                            {
                                const wordList fieldsList(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("sourceFields")); // list of fields to create
                                fvMesh& baffleLessMesh = const_cast<fvMesh&>(meshes_[regioni].time().lookupObject<fvMesh>(meshes_[regioni].name()+".baffleLess"));

                                forAll(fieldsList, fieldi)
                                {
                                   createBaffleLessField<scalar>(fieldsList[fieldi], regioni, baffleLessMesh, nCouplingFields);
                                   createBaffleLessField<vector>(fieldsList[fieldi], regioni, baffleLessMesh, nCouplingFields);
                                   createBaffleLessField<tensor>(fieldsList[fieldi], regioni, baffleLessMesh, nCouplingFields);
                                   createBaffleLessField<symmTensor>(fieldsList[fieldi], regioni, baffleLessMesh, nCouplingFields);
                                   createBaffleLessField<sphericalTensor>(fieldsList[fieldi], regioni, baffleLessMesh, nCouplingFields);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Info<< nl;
}

Foam::dynamicFvMesh& Foam::meshHandler::returnMesh(word meshName)
{
    if (meshName == "dummy")
    {
        return meshes_[0];
    }
    else
    {
        bool found(false);
        label meshIndex(0);
        forAll(meshes_, meshI)
        {
            if (meshes_[meshI].name() == meshName)
            {
                meshIndex = meshI;
                found = true;
                break;
            }
        }

        if (found)
        {
            return meshes_[meshIndex];
        }
        else
        {
            FatalErrorInFunction
                << "Mesh " << meshName << " not found" << nl
                << exit(FatalError);

            return meshes_[0]; //return dummy mesh for compiler
        }
    }
}


Foam::fvMesh& Foam::meshHandler::returnMappingMesh(word meshName)
{
    bool found(false);
    label meshIndex(0);
    forAll(mappingMeshes_,meshI)
    {
        if (mappingMeshes_[meshI].name() == (meshName+".baffleLess"))
        {
            found = true;
            meshIndex = meshI;
            break;
        }
    }

    if (found)
    {
        return mappingMeshes_[meshIndex];
    }
    else
    {
        FatalErrorInFunction
            << "Mesh " << meshName << " not found" << nl
            << exit(FatalError);

        return mappingMeshes_[0]; // return dummymesh for compiler
    }
}



void Foam::meshHandler::interpolateAndMapFields(const Time& runTime)
{

    const IOdictionary couplingDict
    (
        IOobject
        (
            "regionsDict",
            runTime.time().system(),
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if(couplingDict.found("interpolatedMappings"))
    {
        const dictionary interpolatedMappingDict(couplingDict.subDict("interpolatedMappings"));

        forAll(meshes_, regioni)
        {
            if (interpolatedMappingDict.found(meshes_[regioni].name()))
            {
                const dictionary regionToDict(interpolatedMappingDict.subDict(meshes_[regioni].name()));

                scalarList fieldValues(0);
                scalarList xPos(0);
                scalarList yPos(0);
                scalarList zPos(0);

                List<word> regionsFrom(regionToDict.get<List<word>>("sourceRegions"));

                word fieldFromName(regionToDict.get<word>("sourceField"));
                word fieldToName(regionToDict.get<word>("targetField"));

                vector axialDir(regionToDict.get<vector>("axialDirection"));
                scalar magAd = mag(axialDir);
                if (magAd <= VSMALL)
                {
                    FatalErrorInFunction << "axialDirection must be non-zero." << exit(FatalError);
                }
                axialDir /= magAd; // unit axial direction

                scalarList axialLocs(regionToDict.get<List<scalar>>("axialLocations"));

                word interpolationType(regionToDict.get<word>("interpolationType"));

                volScalarField& fieldToBeMapped(const_cast<volScalarField&>(meshes_[regioni].lookupObject<volScalarField>(fieldToName)));


                forAll(regionsFrom, regionFromi)
                {
                    label whichMesh(0);

                    // Get mesh reference
                    forAll(meshes_, i)
                    {
                        if(meshes_[i].name()==regionsFrom[regionFromi])
                            whichMesh = i;
                    }

                    dictionary avgOpts;

                    if (regionToDict.found("averageOptions"))
                    {
                        avgOpts = regionToDict.subDict("averageOptions");
                    }
                    
                    word avgType = avgOpts.get<word>("type");
                    if (avgType != "volumeAverage" and avgType != "patchAverage")
                    {
                        FatalErrorInFunction
                            << "Unknown average type " << avgType
                            << ". Available types: volumeAverage, patchAverage"
                            << exit(FatalError);
                    }

                    scalar dzHalf = 1e-8;
                    if (axialLocs.size() > 1)
                    {
                        scalar minDiff = GREAT;
                        for (auto i = 0; i < axialLocs.size() - 1; ++i)
                        {
                            minDiff = min(minDiff, mag(axialLocs[i+1] - axialLocs[i]));
                        }
                        dzHalf = 0.5*minDiff;
                    }

                    vector centerOfMass
                    (
                        gSum(meshes_[whichMesh].C().field()*meshes_[whichMesh].V().field())
                        / gSum(meshes_[whichMesh].V().field())
                    );

                    const vectorField& C = meshes_[whichMesh].C();
                    const scalarField& V = meshes_[whichMesh].V();
                
                    label patchID = -1;
                    if (avgType == "patchAverage")
                    {
                        word avgPatchName = avgOpts.get<word>("patchName"); // used if patchAverage

                        if (meshes_[whichMesh].boundaryMesh().findPatchID(avgPatchName)==-1)
                        {
                            FatalErrorInFunction << "Patch '" << avgPatchName << "' not found in mesh " << meshes_[whichMesh].name() << exit(FatalError);
                        }
                        patchID = meshes_[whichMesh].boundaryMesh().findPatchID(avgPatchName);
                    }
                
                
                    // --- per-axial-location loop
                    forAll(axialLocs, ai)
                    {
                        scalar sTarget = axialLocs[ai];                    // this is the axial coordinate along axialDir
                        scalar com_s = (centerOfMass & axialDir);         // projection of COM along axialDir

                        // Build a 3D sample point that has projection sTarget on axialDir
                        // i.e. translate the COM along axialDir by (sTarget - com_s)
                        vector sampleVec = centerOfMass + axialDir*(sTarget - com_s);
                        point samplePoint3D(sampleVec);

                        // append global coordinates (so RBF always gets global x,y,z)
                        xPos.append(samplePoint3D.x());
                        yPos.append(samplePoint3D.y());
                        zPos.append(samplePoint3D.z());

                        // compute average value in this axial slice
                        scalar num = 0.0;
                        scalar den = 0.0;
                        scalar avgVal = 0.0;

                        if (avgType == "volumeAverage")
                        {
                            // average over cells whose projection falls inside the slice
                            const volScalarField& srcField = meshes_[whichMesh].lookupObject<volScalarField>(fieldFromName);

                            forAll(C, cellI)
                            {
                                scalar s = (C[cellI] & axialDir);
                                if (mag(s - sTarget) <= dzHalf)
                                {
                                    num += srcField[cellI] * V[cellI];
                                    den += V[cellI];
                                }
                            }

                            // fallback: if slice empty, fallback to nearest cell via findCell(samplePoint3D)
                            if (den <= VSMALL)
                            {
                                label celli = meshes_[whichMesh].findCell(samplePoint3D);
                                if (celli >= 0)
                                {
                                    num = srcField[celli];
                                    den = 1.0;
                                }
                            }
                        }
                        else if (avgType == "patchAverage")
                        {
                            // get patch geometry
                            const polyPatch& pp = meshes_[whichMesh].boundaryMesh()[patchID];
                            const vectorField& faceCentres = pp.faceCentres(); // face centres in global coords
                            const fvPatch& patch = meshes_[whichMesh].boundary()[patchID];
                            const scalarField& faceAreas   = patch.magSf(); // face areas        // face areas

                            scalarField patchVals; 

                            const volScalarField& vField = meshes_[whichMesh].lookupObject<volScalarField>(fieldFromName);
                            patchVals = vField.boundaryField()[patchID];

                            // accumulate area-weighted average of faces inside slice
                            forAll(patchVals, fI)
                            {
                                scalar s = (faceCentres[fI] & axialDir);
                                if (mag(s - sTarget) <= dzHalf)
                                {
                                    num += patchVals[fI] * faceAreas[fI];
                                    den += faceAreas[fI];
                                }
                            }

                            // fallback: if den == 0, try to pick the face closest to samplePoint3D
                            if (den <= VSMALL)
                            {
                                // find nearest face in this patch (cheap linear search)
                                scalar bestDist = GREAT;
                                label bestFace = -1;
                                forAll(faceCentres, fI)
                                {
                                    scalar d = mag(faceCentres[fI] - samplePoint3D);
                                    if (d < bestDist) { bestDist = d; bestFace = fI; }
                                }
                                if (bestFace >= 0)
                                {
                                    num = patchVals[bestFace] * faceAreas[bestFace];
                                    den = faceAreas[bestFace];
                                }
                            }
                        }
                        else
                        {
                            FatalErrorInFunction << "Unknown averageOptions.type: " << avgType << exit(FatalError);
                        }

                        if (den > VSMALL) avgVal = num/den;
                        else avgVal = 0.0; // or handle differently

                        fieldValues.append(avgVal);
                    }

                }
                                    
                scalarList interpolationWeights(0);
                scalar interpolatedValue(0);
                 

                if (interpolationType == "polyharmonicSpline")
                {
                    interpolationWeights = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
                    (
                        xPos, yPos, zPos, fieldValues,  invRBFmatrix_
                    );

                    forAll(meshes_[regioni].C(), centerI)
                    {
                        interpolatedValue = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
                        (
                            interpolationWeights,
                            xPos, yPos, zPos,
                            meshes_[regioni].C()[centerI][0], meshes_[regioni].C()[centerI][1],meshes_[regioni].C()[centerI][2]
                        );

                        fieldToBeMapped[centerI] = interpolatedValue;
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << interpolationType << " is an incorrect "
                        << "radial basis function method. Available methods: polyharmonicSpline, gaussian, kriging"
                        << exit(FatalError);
                }
            }
        }
    }
}


bool Foam::meshHandler::contains(wordList list, word thisWord)
{
    bool found(false);
    forAll(list, wordi)
    {
        if (list[wordi] == thisWord)
        {
            found = true;
        }
    }

    return(found);
}


// ************************************************************************* //
