const word fluidRegionName("fluidRegion");
Foam::fvMesh fluidMesh
(
    Foam::IOobject
    (
        fluidRegionName,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    )
);

autoPtr<fvMesh> fluidMeshNBPtr;

if(runTime.controlDict().lookupOrDefault<bool>("removeBaffles", false))
{
    if (UPstream::nProcs() > 1)
    {
        Info<< "WARNING: the removeBaffles feature is guaranteed to work only "
            << "if each processor contains complete pairs of master an slave "
            << "faces of baffe patches! Otherwise, projection artifacts can "
            << " still occur" << nl << endl;
    }
    
    //- Copy internal geometric data from fluidMesh
    pointField points(fluidMesh.points());
    faceList faces(fluidMesh.faces());
    cellList cells(fluidMesh.cells());

    //- Copy boundary geometric data from fluidMesh
    const polyBoundaryMesh& patches = fluidMesh.boundaryMesh();
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
    fluidMeshNBPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                "fluidRegionNB",
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
    fluidMeshNBPtr().addFvPatches(pList);

    //- Remove baffles
    removeBaffles(fluidMeshNBPtr(), runTime);   
}
else
{
    fluidMeshNBPtr.reset
    (
        new Foam::fvMesh
        (
            Foam::IOobject
            (
                "fluidRegion",
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        )
    );  
}

fvMesh& fluidMeshNB = fluidMeshNBPtr();


//- Silence "unused variable" warnings, I actually do use it and have no idea
//  why the compiler does not "see" it
(void) fluidMeshNB;
