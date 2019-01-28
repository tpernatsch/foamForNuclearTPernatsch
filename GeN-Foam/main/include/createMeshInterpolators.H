
pointField neutroStartPoints = neutroMesh.points(); //used to bring back mesh after interpolations

pointIOField neutroUndeformedPoints
    (
        IOobject
        (
            "neutroUndeformedPoints",
            neutroMesh.time().timeName(),
            neutroMesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        neutroMesh.points()
    );


if(runTime.timeOutputValue()!=0  && fvc::domainIntegrate(neutronics.Disp()).value() != Foam::Vector<double>(0,0,0))
{

    Info << "Bring back meshes to undeformed case to create correct interpolators" << nl << endl; 

    neutroMesh.movePoints(neutroUndeformedPoints);

    Info << "Done" << nl << endl;

}


meshToMesh TMToNeutro
(
    mechMesh,
    neutroMesh,
    Foam::meshToMesh::interpolationMethod(2),
    Foam::meshToMesh::procMapMethod::pmAABB,
    false
);

meshToMesh TMToFluid
(
    mechMesh,
    fluidMeshNB,
    Foam::meshToMesh::interpolationMethod(2),
    Foam::meshToMesh::procMapMethod::pmAABB,
    false
);

meshToMesh neutroToFluid
(
    neutroMesh,
    fluidMeshNB,
    Foam::meshToMesh::interpolationMethod(2),
    Foam::meshToMesh::procMapMethod::pmAABB,
    false
);

meshToMesh fluidToNeutro
(
    fluidMesh,
    neutroMesh,
    Foam::meshToMesh::interpolationMethod(2),
    Foam::meshToMesh::procMapMethod::pmAABB,
    false
);

if(runTime.timeOutputValue()!=0 && fvc::domainIntegrate(neutronics.Disp()).value() != Foam::Vector<double>(0,0,0))
{
    Info << "Re-deform mesh to start calculation" << nl << endl;

    neutroMesh.movePoints(neutroStartPoints);

    neutroMesh.moving(false);

    Info << "Done" << nl << endl;
}



