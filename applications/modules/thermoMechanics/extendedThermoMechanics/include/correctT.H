volScalarField& T(mesh_.lookupObjectRef<volScalarField>("T"));

// This is a post process to match the current structure of GeN-Foam thermomechanics
// The T equation is solved only in the zones where there is no solution from TH
// Here I just set T equal to TStructFromTH if that is the case

forAll(TStructFromTH_, celli)
{
    if(TStructFromTH_[celli] > SMALL)
    {
        T[celli] = TStructFromTH_[celli];
    }
}

