#include "createMeshes.H"

/*
This will soon be used for a multiregion expansion
if (!fluidRegions.size())
{
    FatalErrorIn(args.executable())
        << "No fluid meshes present" << exit(FatalError);
}

fvMesh& mesh = fluidRegions[0];
*/
fvMesh& mesh = fluidMesh;
