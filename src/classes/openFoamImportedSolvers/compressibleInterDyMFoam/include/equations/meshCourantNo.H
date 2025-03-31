    scalar meshCoNum = 0.0;
    scalar meanMeshCoNum = 0.0;

    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(dynMesh.phi()))().primitiveField()
        );

        meshCoNum = 0.5*gMax(sumPhi/dynMesh.V().field())*runTime.deltaTValue();

        meanMeshCoNum =
            0.5*(gSum(sumPhi)/gSum(dynMesh.V().field()))*runTime.deltaTValue();
    }

    Info<< "dynMesh Courant Number mean: " << meanMeshCoNum
    << " max: " << meshCoNum << endl;
 
 // ************************************************************************* //