    

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
if(runTime.controlDict().lookupOrDefault<bool>("removeFluidBaffles", 0))
{
	fluidMeshNBPtr.reset
		(
			new Foam::fvMesh
			(
			    Foam::IOobject
			    (
			        "fluidRegionNB",
			        runTime.timeName(),
			        runTime,
			        Foam::IOobject::MUST_READ
			    )
			)
		);
	removeBaffles(fluidMeshNBPtr(),runTime);	
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


/*
pointField pointsNB = fluidMesh.points();
faceList facesNB = fluidMesh.faces();
labelList ownersNB = fluidMesh.faceOwner();
labelList neighboursNB = fluidMesh.faceOwner();
Foam::fvMesh fluidMeshNB
(
    Foam::IOobject
    (
        "fluidRegionNB",
        runTime.timeName(),
        runTime,
        Foam::IOobject::NO_READ,
        Foam::IOobject::NO_WRITE
    ),
    std::move(pointsNB),
    std::move(facesNB),
    std::move(ownersNB),
    std::move(neighboursNB)
);
*/
//fluidMeshNB.addFvPatches(fluidMesh.);
//removeBaffles(fluidMeshNB,runTime);
