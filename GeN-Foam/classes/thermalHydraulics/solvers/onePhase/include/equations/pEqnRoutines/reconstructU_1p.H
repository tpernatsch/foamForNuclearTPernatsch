//- The ifs are just to avoid useless field reconstructions
//  for integer values of the porousInterfaceSharpness_ factor
if (porousInterfaceSharpness_ != 0.0)
{
    if (porousInterfaceSharpness_ == 1.0)
        U = fvc::reconstruct(MRF_.absolute(fluid_.phi()*alphaf))/alpha;
    else
        U = 
            porousInterfaceSharpness_*
            (
                fvc::reconstruct
                (
                    MRF_.absolute(fluid_.phi()*alphaf)
                )/
                alpha
            )
        +   (1-porousInterfaceSharpness_)*
            (   
                fvc::reconstruct(MRF_.absolute(fluid_.phi()))
            );
}
else
    U = fvc::reconstruct(MRF_.absolute(fluid_.phi()));