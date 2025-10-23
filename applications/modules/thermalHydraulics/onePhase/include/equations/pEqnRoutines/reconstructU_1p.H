//- The ifs are just to avoid useless field reconstructions
//  for integer values of the porousInterfaceSharpness_ factor
if (porousInterfaceSharpness_ != 0.0)
{
    if (porousInterfaceSharpness_ == 1.0)
        U = fvc::reconstruct(phi_)/alpha;
    else
        U = 
            porousInterfaceSharpness_*
            (
                fvc::reconstruct(phi_)/alpha
            )
        +   (1-porousInterfaceSharpness_)*
            (   
                fvc::reconstruct(fluid_.phi())
            );
}
else
{
    // instead of using the fvc::reconstruct, we re-write the reconstruction
    // routine so that the extra patch effects can be corrected when t-jucntion
    // BCs are employed

    surfaceVectorField SfHat(mesh_.Sf()/mesh_.magSf());

    // get the terms of inv(surfaceSum(SfHat*mesh.Sf()))&surfaceSum(SfHat*ssf),

    volTensorField LHS = fvc::surfaceSum(SfHat*mesh_.Sf());

    volVectorField RHS = fvc::surfaceSum(SfHat*fluid_.phi());

    Info << "LHS is " << inv(LHS) << endl;

    U = inv(LHS)&RHS;

}

fluid_.correctAlphaRhoMagU();