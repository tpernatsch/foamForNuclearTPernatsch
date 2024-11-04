Info<< "pTarget_ = " << pTarget_ << nl
    << nl
    << "keff_ previous = " << keff_ << nl
    << endl;

// Compute keff
keff_ *= pTot / pTarget_;

// Limit keff to [0.1, 5]
if (keff_ > 5.0)
{
    Info<< "Warning: keff = " << keff_ << " > 5. Set keff to 5" << endl;
    keff_ = 5.0;
}
else if (keff_ < 0.1)
{
    Info<< "Warning: keff = " << keff_ << " < 0.1. Set keff to 0.1" << endl;
    keff_ = 0.1;
}

Info<< "keff_ = " << keff_ << nl << endl;

// Since a keff is a scalar, the value inside the dictionary is not
// updated at runTime automatically. The line below resets the SCALAR IN
// THE DICTIONARY at runTime so that the dictionary is written with the
// updated value.
reactorState_.set("keff", keff_);
