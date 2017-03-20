Info << "Displace meshes" << nl << endl;

{

	const volPointInterpolation& neutroMeshPointInterpolation = volPointInterpolation::New(neutroMesh);

	pointVectorField neutroPointsDisplacementOld = neutroMeshPointInterpolation.interpolate(DNeutro);

	DNeutro*=0;

	neutroToTM.mapTgtToSrc( meshDisp , plusEqOp<vector>(), DNeutro.primitiveFieldRef());

	DNeutro.correctBoundaryConditions();

	pointVectorField neutroPointsDisplacement = neutroMeshPointInterpolation.interpolate(DNeutro);

	pointField displacedPoints = neutroMesh.points() + neutroPointsDisplacement.internalField() - neutroPointsDisplacementOld.internalField() ;

	neutroMesh.movePoints(displacedPoints);


}
/*
{
	const volPointInterpolation& fluidMeshPointInterpolation = volPointInterpolation::New(fluidMesh);

	pointVectorField fluidPointsDisplacementOld = fluidMeshPointInterpolation.interpolate(DFluid);

	DFluid*=0;

	fluidToTM.mapTgtToSrc( meshDisp , plusEqOp<vector>(), DFluid.primitiveFieldRef());

	DFluid.correctBoundaryConditions();

	pointVectorField fluidPointsDisplacement = fluidMeshPointInterpolation.interpolate(DFluid);

	pointField displacedPoints = fluidMesh.points() + fluidPointsDisplacement.internalField() - fluidPointsDisplacementOld.internalField() ;

	fluidMesh.movePoints(displacedPoints);
}

{
	const volPointInterpolation& mechMeshPointInterpolation = volPointInterpolation::New(mechMesh);

	pointVectorField mechPointsDisplacementOld = mechMeshPointInterpolation.interpolate(meshDispOld);

	meshDisp.correctBoundaryConditions();

	pointVectorField mechPointsDisplacement = mechMeshPointInterpolation.interpolate(meshDisp);

	pointField displacedPoints = mechMesh.points() + mechPointsDisplacement.internalField() - mechPointsDisplacementOld.internalField() ;

	mechMesh.movePoints(displacedPoints);
}
*/



