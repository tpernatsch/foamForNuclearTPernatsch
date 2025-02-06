    scalar meshCoNum = 0.0;
    scalar meanMeshCoNum = 0.0;

    {
        scalarField sumPhi
        (
            fvc::surfaceSum(mag(mesh.phi()))().primitiveField()
        );

        meshCoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();

        meanMeshCoNum =
            0.5*(gSum(sumPhi)/gSum(mesh.V().field()))*runTime.deltaTValue();
    }

    Info<< "mesh Courant Number mean: " << meanMeshCoNum
    << " max: " << meshCoNum << endl;
 
 // ************************************************************************* //