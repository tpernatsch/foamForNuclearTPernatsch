/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2006                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "20/fluidRegion";
    object      T.HX;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 1 0 0 0];


internalField   nonuniform List<scalar> 
120
(
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
1005.04
983.447
966.296
952.851
942.406
934.344
928.152
923.412
919.793
917.035
914.936
913.341
912.128
911.208
910.509
909.978
909.576
909.27
909.038
908.861
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
0
791.839
821.86
843.851
860.076
872.115
881.091
887.807
892.847
896.637
899.493
901.646
903.272
904.501
905.429
906.131
906.662
907.063
907.367
907.596
907.769
0
0
0
0
0
0
0
0
0
0
)
;

boundaryField
{
    outletP
    {
        type            zeroGradient;
    }
    outletS
    {
        type            zeroGradient;
    }
    walls
    {
        type            zeroGradient;
    }
    inletP
    {
        type            zeroGradient;
    }
    inletS
    {
        type            zeroGradient;
    }
}


// ************************************************************************* //
