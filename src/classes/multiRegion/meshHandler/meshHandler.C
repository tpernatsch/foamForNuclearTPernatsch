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

#include "meshHandler.H"
#include "solver.H"
#include "Time.H"
#include "interpolationCellPoint.H"
#include "radialBasisFunctionInterpolation.H"
#include "mergeOrSplitBaffles.H"
#include "hexCellFvMesh.H"



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshHandler::meshHandler(const Time& runTime)
{
    // Read meshHandler dict. Will be used later
    IOdictionary multiRegionCouplingDict
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

    wordList allRegions(controlDict.subDict("regionSolvers").toc());
    wordList meshNames(0);

    forAll(allRegions, regioni)
    {
        word solverName(controlDict.subDict("regionSolvers").get<word>(allRegions[regioni]));
        if (solverName != "multiPhysicsSolver")
        {
            meshNames.append(allRegions[regioni]);
        }
    }

    // If multiPhysicsSolvers are present, also add those extra solvers

    if (multiRegionCouplingDict.found("multiPhysicsSolvers"))
    {
        forAll(multiRegionCouplingDict.subDict("multiPhysicsSolvers").toc(), MPsolvI) // loop over all multiphysics solvers
        {
            wordList subSolvers // get list of solvers belonging to multiphysics solver i
            (
                multiRegionCouplingDict
                    .subDict("multiPhysicsSolvers")
                    .subDict(multiRegionCouplingDict.subDict("multiPhysicsSolvers").toc()[MPsolvI])
                    .subDict("solvers").toc()
            );

            forAll(subSolvers, solvI) // loop ovcer subsolvers and add them to list of meshes to create only if they are not sub-multiphysicssolvers
            {
                if (multiRegionCouplingDict
                    .subDict("multiPhysicsSolvers")
                    .subDict(multiRegionCouplingDict.subDict("multiPhysicsSolvers").toc()[MPsolvI])
                    .subDict("solvers")
                    .get<word>(subSolvers[solvI]) != "multiPhysicsSolver"
                )
                {
                    meshNames.append(subSolvers[solvI]);
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
            new fvMesh
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


            //- Copy internal geometric data from fluidMesh
            pointField points(meshes_[regioni].points());
            faceList faces(meshes_[regioni].faces());
            cellList cells(meshes_[regioni].cells());

            //- Copy boundary geometric data from fluidMesh
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

            //- Assemble mesh as copy from points, faces, cells with no boundary
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

            //- Add boundary data
            baffleLessMesh().addFvPatches(pList);

            //- One of the main issues when passing the actual runTime to removeBaffles
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

            //- Copy mesh zones
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
            "multiRegionCouplingDict",
            runTime.time().constant(),
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
                    if (regionsInDict[i]!=meshes_[regioni].name()) //-look for other (than regionMehses_[regioni]) regions source fields
                    {
                        const dictionary regionFromDict(mappingDict.subDict(regionsInDict[i]));
                        const wordList regionsFrom(regionFromDict.toc());
                        // scalarCouplingFields_[regioni].setSize(regionsFrom.size());

                        forAll(regionsFrom, regionFromi)//-loop over regions in the sub-dict
                        {
                            if(regionsFrom[regionFromi] == meshes_[regioni].name())//-if i find regionMehses_[regioni] in the subdict
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

Foam::fvMesh& Foam::meshHandler::returnMesh(word meshName)
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



// void Foam::meshHandler::interpolateAndMapFields(const Time& runTime)
// {

//     const IOdictionary couplingDict
//     (
//         IOobject
//         (
//             "multiRegionCouplingDict",
//             runTime.time().constant(),
//             runTime.db(),
//             IOobject::MUST_READ,
//             IOobject::NO_WRITE
//         )
//     );

//     if(couplingDict.found("interpolatedMappings"))
//     {
//         const dictionary interpolatedMappingDict(couplingDict.subDict("interpolatedMappings"));

//         forAll(meshes_, regioni)
//         {
//             if (interpolatedMappingDict.found(meshes_[regioni].name()))
//             {
//                 const dictionary regionToDict(interpolatedMappingDict.subDict(meshes_[regioni].name()));

//                 scalarList fieldValues(0);
//                 scalarList xPos(0);
//                 scalarList yPos(0);
//                 scalarList zPos(0);

//                 List<word> regionsFrom(regionToDict.get<List<word>>("fromWhichRegions"));
//                 word fieldFromName(regionToDict.get<word>("fieldFromName"));
//                 word fieldToName(regionToDict.get<word>("fieldToName"));
//                 scalarList axialLocs(regionToDict.get<List<scalar>>("axialLocations"));
//                 word interpolationType(regionToDict.get<word>("interpolationType"));
//                 volScalarField& fieldToBeMapped(const_cast<volScalarField&>(meshes_[regioni].lookupObject<volScalarField>(fieldToName)));


//                 forAll(regionsFrom, regionFromi)
//                 {
//                     label whichMesh(0);

//                     // Get mesh reference
//                     forAll(meshes_, i)
//                     {
//                         if(meshes_[i].name()==regionsFrom[regionFromi])
//                             whichMesh = i;
//                     }


//                     // - Get axial locations
//                     zPos.append(axialLocs);

//                     // - Get x and y
//                     vector centerOfMass(gSum(meshes_[whichMesh].C().field()*meshes_[whichMesh].V().field())/gSum(meshes_[whichMesh].V().field()));

//                     for(label i = 0; i<axialLocs.size(); i++)
//                     {
//                         xPos.append(centerOfMass[0]);
//                         yPos.append(centerOfMass[1]);
//                         // - Now I have a list of (x,y,z) for one region. I need to get to associate a value to each coordinate

//                         point samplePoint(centerOfMass[0], centerOfMass[1], axialLocs[i]);

//                         // interpolationCellPoint<scalar> pointInterpolator(meshes_[whichMesh]);

//                         label celli = meshes_[whichMesh].findCell(samplePoint);

//                         const volScalarField& field(meshes_[whichMesh].lookupObject<volScalarField>(fieldFromName));

//                         fieldValues.append(field[celli]);

//                     }
//                 }

//                 scalarList interpolationWeights(0);
//                 scalar interpolatedValue(0);


//                 if (interpolationType == "kriging")
//                 {
//                     Foam::radialBasisFunctionInterpolation::solveKriging
//                     (
//                         xPos, yPos, zPos, invRBFmatrix_
//                     );

//                     forAll(meshes_[regioni].C(), centerI)
//                     {
//                         interpolatedValue = Foam::radialBasisFunctionInterpolation::kriging
//                         (
//                             xPos, yPos, zPos, fieldValues,
//                             meshes_[regioni].C()[centerI][0], meshes_[regioni].C()[centerI][1],meshes_[regioni].C()[centerI][2],
//                             invRBFmatrix_
//                         );

//                         fieldToBeMapped[centerI] = interpolatedValue;
//                     }

//                     fieldToBeMapped.correctBoundaryConditions();
//                 }
//                 else if (interpolationType == "polyharmonicSpline")
//                 {
//                     interpolationWeights = Foam::radialBasisFunctionInterpolation::solvePolyharmonicSpline
//                     (
//                         xPos, yPos, zPos, fieldValues,  invRBFmatrix_
//                     );

//                     forAll(meshes_[regioni].C(), centerI)
//                     {
//                         interpolatedValue = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
//                         (
//                             interpolationWeights,
//                             xPos, yPos, zPos,
//                             meshes_[regioni].C()[centerI][0], meshes_[regioni].C()[centerI][1],meshes_[regioni].C()[centerI][2]
//                         );

//                         fieldToBeMapped[centerI] = interpolatedValue;
//                     }
//                 }
//                 else
//                 {
//                     FatalErrorInFunction
//                         << interpolationType << " is an incorrect "
//                         << "radial basis function method. Available methods: polyharmonicSpline, gaussian, kriging"
//                         << exit(FatalError);
//                 }
//             }
//         }
//     }
// }


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
