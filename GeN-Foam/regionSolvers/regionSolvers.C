/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
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

#include "regionSolvers.H"
#include "solver.H"
#include "Time.H"
#include "interpolationCellPoint.H"
#include "radialBasisFunctionInterpolation.H"
#include "mergeOrSplitBaffles.H"
#include "fluid.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSolvers::regionSolvers(const Time& runTime)
{
    List<Pair<word>> regionSolverNames;


    //Read mapping dict. Will be used later
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


    if (controlDict.found("regionSolvers"))
    {
        const dictionary& regionSolversDict =
            controlDict.subDict("regionSolvers");

        forAllConstIter(dictionary, regionSolversDict, iter)
        {
            const word regionName(iter().keyword());
            const word solverName(iter().stream());

            regionSolverNames.append(Pair<word>(regionName, solverName));
        }
    }
    else //no backward compatibility for now
    {
        FatalIOErrorInFunction(runTime.controlDict())
                        << "regionSolvers list missing from "
                        << runTime.controlDict().name()
                        << exit(FatalIOError);
    }


    regionMeshes_.setSize(regionSolverNames.size());

    //-Set size of list of (possibly) baffleLess meshes
    mappingRegionMeshes_.setSize(regionSolverNames.size());
    //-Set size of list of couplingFields
    scalarCouplingFields_.setSize(regionSolverNames.size());
    vectorCouplingFields_.setSize(regionSolverNames.size());

    solvers_.setSize(regionSolverNames.size());
    prefixes_.setSize(regionSolverNames.size());

    string::size_type nRegionNameChars = 0;

    forAll(regionSolverNames, i)
    {
        const word& regionName = regionSolverNames[i].first();
        const word& solverName = regionSolverNames[i].second();

        // Load the solver library -> not sure it is necessary for now
        // solver::load(solverName);
        regionMeshes_.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.name(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );

        Info << "Creating solver for region "<< regionName<<nl<<endl;
        solvers_.set(i, solver::New(solverName, regionMeshes_[i]));

        prefixes_[i] = regionName;
        nRegionNameChars = max(nRegionNameChars, regionName.size());
    }


    createBaffleLessMeshes(runTime);
    createCouplingFields(runTime);

    nRegionNameChars++;

    prefix0_.append(nRegionNameChars, ' ');

    forAll(regionSolverNames, i)
    {
        prefixes_[i].append(nRegionNameChars - prefixes_[i].size(), ' ');
    }


    //Create mappings at beginning


    multiRegionResidual_=(multiRegionCouplingDict.get<scalar>("multiRegionResidual"));
    maxMultiRegionIterations_=(multiRegionCouplingDict.get<scalar>("maxMultiRegionIterations"));
    multiRegionTightCoupled_=(multiRegionCouplingDict.getOrDefault<bool>("multiRegionTightCoupled", false));
    // Get info about fields mapping
    mappingDict_=multiRegionCouplingDict.subDict("mappings");


    mappingList_.setSize(regionMeshes_.size());

    forAll(regionMeshes_,i)
    {
        mappingList_.set
        (
            i,
            new PtrList<meshToMesh>
            (
                regionMeshes_.size()
            )
        );
    }

    //Loop over all regions
    forAll(regionMeshes_, i)
    {
        if (mappingDict_.found(regionMeshes_[i].name())) // check if i find entry in the mapping dict
        {
            forAll(regionMeshes_, j) // loop over all other regions
            {
                if (mappingDict_.subDict(regionMeshes_[i].name()).found(regionMeshes_[j].name())) // check if i find entry in the mapping subDict
                {
                    mappingList_[i].set
                    (
                        j,
                        new meshToMesh
                        (
                            regionMeshes_[i],
                            controlDict.subDict("removeBaffles").get<bool>(regionMeshes_[j].name())?
                            mappingRegionMeshes_[j]
                            :regionMeshes_[j],
                            Foam::meshToMesh::interpolationMethod::imCellVolumeWeight, // for now hard-coded, possibly implement mapping type via dict
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

Foam::regionSolvers::~regionSolvers()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionSolvers::setGlobalPrefix() const
{
    Pout.prefix() = prefix0_;
}


void Foam::regionSolvers::setPrefix(const label i) const
{
    Pout.prefix() = prefixes_[i];
}


void Foam::regionSolvers::resetPrefix() const
{
    Pout.prefix() = string::null;
}

void Foam::regionSolvers::mapFields( const Time& runTime)
{

    forAll(regionMeshes_, i)
    {
        forAll(regionMeshes_, j)
        {
            if(mappingDict_.subDict(regionMeshes_[i].name()).found(regionMeshes_[j].name()))
            {
                dictionary iTojRegionmappingDict_(mappingDict_.subDict(regionMeshes_[i].name()).subDict(regionMeshes_[j].name()));
                List<word> targetFields(iTojRegionmappingDict_.getOrDefault<List<word>>("targetFields", List<word>::null()));
                List<word> sourceFields(iTojRegionmappingDict_.getOrDefault<List<word>>("sourceFields",List<word>::null()));
                List<word> fieldTypes(iTojRegionmappingDict_.getOrDefault<List<word>>("fieldTypes",List<word>::null()));

                bool removeBaffles(false);

                if(runTime.controlDict().found("removeBaffles"))
                {
                    const dictionary& removeBafflesDict(runTime.controlDict().subDict("removeBaffles"));
            
                    if(removeBafflesDict.getOrDefault<bool>(regionMeshes_[j].name(), false))
                    {
                        forAll(sourceFields, sourcei)
                        {
                            sourceFields[sourcei]+=".baffleLess";
                        }
                        removeBaffles = true;
                    }
                }

                if(targetFields.size()!=sourceFields.size())
                {
                    FatalErrorInFunction
                    << "Number of target fields from region "<< regionMeshes_[j].name()<< " is not equal to number of source fields to region "<<mappingRegionMeshes_[i].name()<<"!"
                    << exit(FatalError);
                }
                else if(fieldTypes.size()!=sourceFields.size())
                {
                    FatalErrorInFunction
                    << "Number of specified field types in mapping from region "<< regionMeshes_[j].name()<< " to region "<<regionMeshes_[i].name()<<" is not adequate!"
                    << exit(FatalError);
                }
                else
                {
                    forAll(sourceFields, fieldi)
                    {

                        if(fieldTypes[fieldi]=="scalar")
                        {
                            volScalarField& tgtField = const_cast<volScalarField&>(regionMeshes_[i].lookupObject<volScalarField>(targetFields[fieldi]));
                            if ((removeBaffles and !(mappingRegionMeshes_[j].foundObject<volScalarField>(sourceFields[fieldi]))) or (!removeBaffles and !(regionMeshes_[j].foundObject<volScalarField>(sourceFields[fieldi]))))
                            {
                                Info<<"Warning! Field " <<sourceFields[fieldi]<< " not found! "<<endl;

                            }
                            else
                            {
                                mappingList_[i][j].mapTgtToSrc
                                (
                                    removeBaffles?
                                    mappingRegionMeshes_[j].lookupObject<volScalarField>(sourceFields[fieldi])
                                    :regionMeshes_[j].lookupObject<volScalarField>(sourceFields[fieldi]),
                                    plusEqOp<scalar>(),
                                    tgtField
                                );
                                tgtField.correctBoundaryConditions();
                            }
                        }
                        else if(fieldTypes[fieldi]=="vector")
                        {
                            volVectorField& tgtField = const_cast<volVectorField&>(regionMeshes_[i].lookupObject<volVectorField>(targetFields[fieldi]));
                            mappingList_[i][j].mapTgtToSrc
                            (
                                removeBaffles?
                                mappingRegionMeshes_[j].lookupObject<volVectorField>(sourceFields[fieldi])
                                :regionMeshes_[j].lookupObject<volVectorField>(sourceFields[fieldi]),
                                plusEqOp<vector>(),
                                tgtField
                            );
                            tgtField.correctBoundaryConditions();
                        }
                        else
                        {
                            FatalErrorInFunction
                            << "Field type " << fieldTypes[fieldi]<< " in mapping from region "<< regionMeshes_[j].name()<< " to region "<<mappingRegionMeshes_[i].name()<<" is not known!"
                            << exit(FatalError);
                        }
                    }
                }
            }
        }
    }

}


void Foam::regionSolvers::initializeMappedFields( const Time& runTime)
{

    forAll(regionMeshes_, i)
    {
        forAll(regionMeshes_, j)
        {
            if(mappingDict_.subDict(regionMeshes_[i].name()).found(regionMeshes_[j].name()))
            {
                dictionary iTojRegionmappingDict_(mappingDict_.subDict(regionMeshes_[i].name()).subDict(regionMeshes_[j].name()));
                List<word> targetFields(iTojRegionmappingDict_.getOrDefault<List<word>>("targetFields", List<word>::null()));
                List<word> sourceFields(iTojRegionmappingDict_.getOrDefault<List<word>>("sourceFields",List<word>::null()));
                List<word> fieldTypes(iTojRegionmappingDict_.getOrDefault<List<word>>("fieldTypes",List<word>::null()));

                bool removeBaffles(false);

                if(runTime.controlDict().found("removeBaffles"))
                {
                    const dictionary& removeBafflesDict(runTime.controlDict().subDict("removeBaffles"));
            
                    if(removeBafflesDict.getOrDefault<bool>(regionMeshes_[j].name(), false))
                    {
                        forAll(sourceFields, sourcei)
                        {
                            sourceFields[sourcei]+=".baffleLess";
                        }
                        removeBaffles = true;
                    }
                }

                if(targetFields.size()!=sourceFields.size())
                {
                    FatalErrorInFunction
                    << "Number of target fields from region "<< regionMeshes_[j].name()<< " is not equal to number of source fields to region "<<mappingRegionMeshes_[i].name()<<"!"
                    << exit(FatalError);
                }
                else if(fieldTypes.size()!=sourceFields.size())
                {
                    FatalErrorInFunction
                    << "Number of specified field types in mapping from region "<< regionMeshes_[j].name()<< " to region "<<regionMeshes_[i].name()<<" is not adequate!"
                    << exit(FatalError);
                }
                else
                {
                    forAll(sourceFields, fieldi)
                    {

                        if(fieldTypes[fieldi]=="scalar")
                        {
                            volScalarField& tgtField = const_cast<volScalarField&>(regionMeshes_[i].lookupObject<volScalarField>(targetFields[fieldi]));
                            if ((removeBaffles and !(mappingRegionMeshes_[j].foundObject<volScalarField>(sourceFields[fieldi]))) or (!removeBaffles and !(regionMeshes_[j].foundObject<volScalarField>(sourceFields[fieldi]))))
                            {
                                Info<<"Warning! Field " <<sourceFields[fieldi]<< " not found! "<<endl;

                            }
                            else
                            {
                                mappingList_[i][j].mapTgtToSrc
                                (
                                    removeBaffles?
                                    mappingRegionMeshes_[j].lookupObject<volScalarField>(sourceFields[fieldi])
                                    :regionMeshes_[j].lookupObject<volScalarField>(sourceFields[fieldi]),
                                    plusEqOp<scalar>(),
                                    tgtField
                                );
                                tgtField.correctBoundaryConditions();
                                tgtField.write();
                            }
                        }
                        else if(fieldTypes[fieldi]=="vector")
                        {
                            volVectorField& tgtField = const_cast<volVectorField&>(regionMeshes_[i].lookupObject<volVectorField>(targetFields[fieldi]));
                            mappingList_[i][j].mapTgtToSrc
                            (
                                removeBaffles?
                                mappingRegionMeshes_[j].lookupObject<volVectorField>(sourceFields[fieldi])
                                :regionMeshes_[j].lookupObject<volVectorField>(sourceFields[fieldi]),
                                plusEqOp<vector>(),
                                tgtField
                            );
                            tgtField.correctBoundaryConditions();
                            tgtField.write();
                        }
                        else
                        {
                            FatalErrorInFunction
                            << "Field type " << fieldTypes[fieldi]<< " in mapping from region "<< regionMeshes_[j].name()<< " to region "<<mappingRegionMeshes_[i].name()<<" is not known!"
                            << exit(FatalError);
                        }
                    }
                }
            }
        }
    }

}



void Foam::regionSolvers::createBaffleLessMeshes(const Time& runTime)
{
   
    forAll(regionMeshes_, regioni)
    {        
        autoPtr<fvMesh> baffleLessMesh;

        bool removeBaffleBool(false);

        if (runTime.controlDict().found("removeBaffles"))
        {
            const dictionary& removeBafflesDict = runTime.controlDict().subDict("removeBaffles");
            removeBaffleBool = removeBafflesDict.getOrDefault<bool>(regionMeshes_[regioni].name(), "false");
        }

        if (removeBaffleBool)
        {

            Info << "Creating baffleless mesh for region " << regionMeshes_[regioni].name()<<endl;
            if (UPstream::nProcs() > 1)
            {
                Info<< "WARNING: the removeBaffles feature is guaranteed to work only "
                    << "if each processor contains complete pairs of master and slave "
                    << "faces of baffe patches! Otherwise, projection artifacts can "
                    << " still occur" << nl << endl;
            }

                
            //- Copy internal geometric data from fluidMesh
            pointField points(regionMeshes_[regioni].points());
            faceList faces(regionMeshes_[regioni].faces());
            cellList cells(regionMeshes_[regioni].cells());

            //- Copy boundary geometric data from fluidMesh
            const polyBoundaryMesh& patches = regionMeshes_[regioni].boundaryMesh();
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
                        IOobject::groupName(regionMeshes_[regioni].name(), "baffleLess"),
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
                new Time(Foam::Time::controlDictName, regionMeshes_[regioni].time().rootPath(), regionMeshes_[regioni].time().caseName())
            );
            removeBaffles(baffleLessMesh(), dummyRunTimePtr());
            dummyRunTimePtr.clear();    

            //- Copy mesh zones
            List<pointZone*> pointZonesNB(0);
            List<faceZone*> faceZonesNB(0);
            List<cellZone*> cellZonesNB(0);
            forAll(regionMeshes_[regioni].pointZones().names(), i)
            {
                word name(regionMeshes_[regioni].pointZones().names()[i]);
                const pointZone& origZone(regionMeshes_[regioni].pointZones()[name]);
                pointZonesNB.append
                (
                    new pointZone
                    (
                        name,
                        origZone,
                        i,
                        regionMeshes_[regioni].pointZones()
                    )
                );
            }
            forAll(regionMeshes_[regioni].faceZones().names(), i)
            {
                word name(regionMeshes_[regioni].faceZones().names()[i]);
                const faceZone& origZone(regionMeshes_[regioni].faceZones()[name]);
                faceZonesNB.append
                (
                    new faceZone
                    (
                        name,
                        origZone,
                        origZone.flipMap(),
                        i,
                        regionMeshes_[regioni].faceZones()
                    )
                );
            }
            forAll(regionMeshes_[regioni].cellZones().names(), i)
            {
                word name(regionMeshes_[regioni].cellZones().names()[i]);
                const cellZone& origZone(regionMeshes_[regioni].cellZones()[name]);
                cellZonesNB.append
                (
                    new cellZone
                    (
                        name,
                        origZone,
                        i,
                        regionMeshes_[regioni].cellZones()
                    )
                );
            }
            baffleLessMesh().addZones(pointZonesNB, faceZonesNB, cellZonesNB);

        }

        else
        {
            baffleLessMesh = nullptr;
        }
        

        mappingRegionMeshes_.set(regioni,baffleLessMesh);
    }

}



void Foam::regionSolvers::createCouplingFields(const Time& runTime)
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

    forAll(regionMeshes_, regioni) //iterate over the meshes
    {       
        if(runTime.controlDict().found("removeBaffles"))
        {
            const dictionary& removeBafflesDict = runTime.controlDict().subDict("removeBaffles");
            if(removeBafflesDict.getOrDefault<bool>(regionMeshes_[regioni].name(),false))
            {

                label nScalarFields(0);
                label nVectorFields(0);
                
                scalarCouplingFields_.set
                (
                    regioni,
                    new PtrList<volScalarField>(nScalarFields+1)
                );
                vectorCouplingFields_.set
                (
                    regioni,
                    new PtrList<volVectorField>(nVectorFields+1)
                );
                
                forAll(regionsInDict, i) //iterate over the regions found in the dictionary
                {
                
                    if(regionsInDict[i]!=regionMeshes_[regioni].name()) //-look for other (than regionMehses_[regioni]) regions source fields
                    {
                        
                        const dictionary regionFromDict(mappingDict.subDict(regionsInDict[i]));
                        const wordList regionsFrom(regionFromDict.toc());
                        // scalarCouplingFields_[regioni].setSize(regionsFrom.size());
                        
                        forAll(regionsFrom, regionFromi)//-loop over regions in the sub-dict
                        {
                            
                            if(regionsFrom[regionFromi]==regionMeshes_[regioni].name())//-if i find regionMehses_[regioni] in the subdict
                            {


                                const wordList fieldsList(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("sourceFields")); // list of fields to create
                                const wordList fieldTypes(regionFromDict.subDict(regionsFrom[regionFromi]).get<wordList>("fieldTypes")); // list of types of fields to create
                                fvMesh& baffleLessMesh = const_cast<fvMesh&>(regionMeshes_[regioni].time().lookupObject<fvMesh>(regionMeshes_[regioni].name()+".baffleLess"));

                                forAll(fieldsList, fieldi)
                                {
                                    if (fieldTypes[fieldi] == "scalar" and (!baffleLessMesh.foundObject<volScalarField>(fieldsList[fieldi]+".baffleLess")))
                                    {
                                        scalarCouplingFields_[regioni].resize(nScalarFields+1);
                                        scalarCouplingFields_[regioni].set
                                        (
                                            nScalarFields,
                                            new volScalarField
                                            (
                                                IOobject
                                                (
                                                    IOobject::groupName(fieldsList[fieldi], "baffleLess"),
                                                    baffleLessMesh.time().name(),
                                                    baffleLessMesh,
                                                    IOobject::NO_READ,
                                                    IOobject::NO_WRITE
                                                ),
                                                baffleLessMesh,
                                                dimensionedScalar("", regionMeshes_[regioni].lookupObject<volScalarField>(fieldsList[fieldi]).dimensions() ,0.0),
                                                zeroGradientFvPatchScalarField::typeName
                                            )
                                        );
                                        nScalarFields++;

                                        Info << "Creating baffleless field "<<scalarCouplingFields_[regioni][nScalarFields-1].name() <<endl;
                                    }
                                    else if (fieldTypes[fieldi] == "vector" and (!baffleLessMesh.foundObject<volVectorField>(fieldsList[fieldi]+".baffleLess")))
                                    {
                                        vectorCouplingFields_[regioni].resize(nVectorFields+1);
                                        vectorCouplingFields_[regioni].set
                                        (
                                            nVectorFields,
                                            new volVectorField
                                            (
                                                IOobject
                                                (
                                                    IOobject::groupName(fieldsList[fieldi], "baffleLess"),
                                                    baffleLessMesh.time().name(),
                                                    baffleLessMesh,
                                                    IOobject::NO_READ,
                                                    IOobject::NO_WRITE
                                                ),
                                                baffleLessMesh,
                                                dimensionedVector("", regionMeshes_[regioni].lookupObject<volVectorField>(fieldsList[fieldi]).dimensions() ,vector::zero),
                                                zeroGradientFvPatchScalarField::typeName
                                            )
                                        );
                                        nVectorFields++;

                                        Info << "Creatinfg baffleless field "<< vectorCouplingFields_[regioni][nVectorFields-1].name() <<endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Info <<nl;
}


// void Foam::regionSolvers::interpolateAndMapFields(const Time& runTime)
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

//         forAll(regionMeshes_, regioni)
//         {
//             if (interpolatedMappingDict.found(regionMeshes_[regioni].name()))
//             {
//                 const dictionary regionToDict(interpolatedMappingDict.subDict(regionMeshes_[regioni].name()));

//                 scalarList fieldValues(0);
//                 scalarList xPos(0);
//                 scalarList yPos(0);
//                 scalarList zPos(0);

//                 List<word> regionsFrom(regionToDict.get<List<word>>("fromWhichRegions"));
//                 word fieldFromName(regionToDict.get<word>("fieldFromName"));
//                 word fieldToName(regionToDict.get<word>("fieldToName"));
//                 scalarList axialLocs(regionToDict.get<List<scalar>>("axialLocations"));
//                 word interpolationType(regionToDict.get<word>("interpolationType"));
//                 volScalarField& fieldToBeMapped(const_cast<volScalarField&>(regionMeshes_[regioni].lookupObject<volScalarField>(fieldToName)));


//                 forAll(regionsFrom, regionFromi)
//                 {
//                     label whichMesh(0);

//                     // Get mesh reference
//                     forAll(regionMeshes_, i)
//                     {
//                         if(regionMeshes_[i].name()==regionsFrom[regionFromi])
//                             whichMesh = i;
//                     }


//                     // - Get axial locations
//                     zPos.append(axialLocs);

//                     // - Get x and y
//                     vector centerOfMass(gSum(regionMeshes_[whichMesh].C().field()*regionMeshes_[whichMesh].V().field())/gSum(regionMeshes_[whichMesh].V().field()));

//                     for(label i = 0; i<axialLocs.size(); i++)
//                     {
//                         xPos.append(centerOfMass[0]);
//                         yPos.append(centerOfMass[1]);
//                         // - Now I have a list of (x,y,z) for one region. I need to get to associate a value to each coordinate

//                         point samplePoint(centerOfMass[0], centerOfMass[1], axialLocs[i]);

//                         // interpolationCellPoint<scalar> pointInterpolator(regionMeshes_[whichMesh]);

//                         label celli = regionMeshes_[whichMesh].findCell(samplePoint);
                        
//                         const volScalarField& field(regionMeshes_[whichMesh].lookupObject<volScalarField>(fieldFromName));

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

//                     forAll(regionMeshes_[regioni].C(), centerI)
//                     {
//                         interpolatedValue = Foam::radialBasisFunctionInterpolation::kriging
//                         (
//                             xPos, yPos, zPos, fieldValues,
//                             regionMeshes_[regioni].C()[centerI][0], regionMeshes_[regioni].C()[centerI][1],regionMeshes_[regioni].C()[centerI][2],
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

//                     forAll(regionMeshes_[regioni].C(), centerI)
//                     {
//                         interpolatedValue = Foam::radialBasisFunctionInterpolation::polyharmonicSpline
//                         (
//                             interpolationWeights,
//                             xPos, yPos, zPos,
//                             regionMeshes_[regioni].C()[centerI][0], regionMeshes_[regioni].C()[centerI][1],regionMeshes_[regioni].C()[centerI][2]
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



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::solver& Foam::regionSolvers::operator[](const label i)
{
    setPrefix(i);
    return solvers_[i];
}


// ************************************************************************* //
