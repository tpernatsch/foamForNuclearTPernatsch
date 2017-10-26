Info << "Displace meshes" << nl << endl;

{
    const volPointInterpolation& neutroMeshPointInterpolation = volPointInterpolation::New(neutroMesh);

    tmp<pointVectorField> neutroPointsDisplacementOld = neutroMeshPointInterpolation.interpolate(DNeutro);

    DNeutro*=0;
    neutroToTM.mapTgtToSrc( meshDisp , plusEqOp<vector>(), DNeutro.primitiveFieldRef());
    DNeutro.correctBoundaryConditions();

    tmp<pointVectorField> neutroPointsDisplacement = neutroMeshPointInterpolation.interpolate(DNeutro);

    tmp<pointField> displacedPoints = neutroMesh.points()
                                    + neutroPointsDisplacement->internalField()
                                    - neutroPointsDisplacementOld->internalField() ;

    neutroMesh.movePoints(displacedPoints);

}

{
    const volPointInterpolation& fluidMeshPointInterpolation = volPointInterpolation::New(fluidMesh);

    tmp<pointVectorField> fluidPointsDisplacementOld = fluidMeshPointInterpolation.interpolate(DFluid);

    DFluid*=0;
    fluidToTM.mapTgtToSrc( meshDisp , plusEqOp<vector>(), DFluid.primitiveFieldRef());
    DFluid.correctBoundaryConditions();

    tmp<pointVectorField> fluidPointsDisplacement = fluidMeshPointInterpolation.interpolate(DFluid);

    tmp<pointField> displacedPoints = fluidMesh.points() 
                                    + fluidPointsDisplacement->internalField()
                                    - fluidPointsDisplacementOld->internalField() ;

    fluidMesh.movePoints(displacedPoints);

} 

{
    const volPointInterpolation& mechMeshPointInterpolation = volPointInterpolation::New(mechMesh);

    tmp<pointVectorField> mechPointsDisplacementOld = mechMeshPointInterpolation.interpolate(meshDisp.oldTime());

    meshDisp.correctBoundaryConditions();

    tmp<pointVectorField> mechPointsDisplacement = mechMeshPointInterpolation.interpolate(meshDisp);

    tmp<pointField> displacedPoints = mechMesh.points()
                                    + mechPointsDisplacement->internalField()
                                    - mechPointsDisplacementOld->internalField() ;

    mechMesh.movePoints(displacedPoints);

}






