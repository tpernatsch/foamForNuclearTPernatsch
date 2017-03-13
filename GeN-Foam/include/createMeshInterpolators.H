
pointField neutroStartPoints = neutroMesh.points(); //used to bring back mesh after interpolations
pointField fluidStartPoints = fluidMesh.points(); //used to bring back mesh after interpolations
pointField mechStartPoints = mechMesh.points(); //used to bring back mesh after interpolations


if(runTime.timeOutputValue()!=0  && fvc::domainIntegrate(DNeutro).value() != Foam::Vector<double>(0,0,0))
{
	Info << "Bring back meshes to underformed case to create correct interpolators" << nl << endl;

	const volPointInterpolation& neutroMeshPointInterpolationInitial = volPointInterpolation::New(neutroMesh);

	pointVectorField neutroPointsDisplacementInitial = neutroMeshPointInterpolationInitial.interpolate(DNeutro);

	pointField neutroDisplacedPointsInitial = neutroStartPoints - neutroPointsDisplacementInitial.internalField() ;

	neutroMesh.movePoints(neutroDisplacedPointsInitial);



	const volPointInterpolation& fluidMeshPointInterpolationInitial = volPointInterpolation::New(fluidMesh);

	pointVectorField fluidPointsDisplacementInitial = fluidMeshPointInterpolationInitial.interpolate(DFluid);

	pointField fluidDisplacedPointsInitial = fluidStartPoints - fluidPointsDisplacementInitial.internalField() ;

	fluidMesh.movePoints(fluidDisplacedPointsInitial);



	const volPointInterpolation& mechMeshPointInterpolationInitial = volPointInterpolation::New(mechMesh);

	pointVectorField mechPointsDisplacementInitial = mechMeshPointInterpolationInitial.interpolate(meshDisp);

	pointField mechDisplacedPointsInitial = mechStartPoints - mechPointsDisplacementInitial.internalField() ;

	mechMesh.movePoints(mechDisplacedPointsInitial);

	Info << "Done" << nl << endl;
}


	meshToMesh neutroToTM
	(
		neutroMesh,
		mechMesh,
		Foam::meshToMesh::interpolationMethod(2),
		false
	);

	meshToMesh fluidToTM
	(
		fluidMesh,
		mechMesh,
		Foam::meshToMesh::interpolationMethod(2),
		false
	);


	meshToMesh neutroToFluid
	(
		neutroMesh,
		fluidMesh,
		Foam::meshToMesh::interpolationMethod(2),
		false
	);


if(runTime.timeOutputValue()!=0 && fvc::domainIntegrate(DNeutro).value() != Foam::Vector<double>(0,0,0))
{
	Info << "Re-derform mesh to start calculation" << nl << endl;

	neutroMesh.movePoints(neutroStartPoints);
	fluidMesh.movePoints(fluidStartPoints);
	mechMesh.movePoints(mechStartPoints);
/*
	neutroMesh.Foam::polyMesh::resetMotion();
	fluidMesh.Foam::polyMesh::resetMotion();
	mechMesh.Foam::polyMesh::resetMotion();
*/

	neutroMesh.readUpdate();
	fluidMesh.readUpdate();
	mechMesh.readUpdate();

	neutroMesh.moving(false);
	fluidMesh.moving(false);
	mechMesh.moving(false);

	Info << "Done" << nl << endl;
}















