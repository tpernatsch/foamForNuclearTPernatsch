// Rescales Keff_ to drive the computed power Q_ toward the target
// Qnominal_, exactly mirroring applications/modules/neutronics's own
// calcKeff.H. Since updateCrossSections.H divides Sf_ by Keff_, this acts
// as a proportional controller: too much power raises Keff_ (suppressing
// the next iteration's fission source), too little lowers it.

{
    scalar pTot = max(Q_).value();
    scalar pTarget = Qnominal_.value();

    Info<< "Target power = " << pTarget/1e6 << " MW" << nl
        << "Keff previous = " << max(Keff_).value() << endl;

    Keff_ *= pTot/pTarget;

    if (max(Keff_).value() > 5.0)
    {
        Info<< "Warning: Keff = " << max(Keff_).value() << " > 5. Set Keff to 5" << endl;
        Keff_ = dimensionedScalar("", dimless, 5.0);
    }
    else if (max(Keff_).value() < 0.1)
    {
        Info<< "Warning: Keff = " << max(Keff_).value() << " < 0.1. Set Keff to 0.1" << endl;
        Keff_ = dimensionedScalar("", dimless, 0.1);
    }

    Info<< "Keff = " << max(Keff_).value() << nl << endl;

    reactivity_ = (Keff_ - 1)/Keff_;
}
