{
    // Reset fields
    qPrompt_ *= 0.0;
    qDecay_ *= 0.0;

    // Evaluate prompt heat source from current fluxes and cross-sections
    forAll(flux_, i)
    {
        qPrompt_ += (1 - decBetaTot_)*Ef_[i]*Sf_[i]*flux_[i];
    }

    // Evaluate delayed heat source from current precursors
    forAll(dec_, i)
    {
        qDecay_ += decLambda_[i]*dec_[i];
    }

    q_ = qPrompt_ + qDecay_;

    Qprompt_ = unitField_*fvc::domainIntegrate(qPrompt_).value()*dimensionedScalar("", dimPower, 1);

    Qdecay_ = unitField_*fvc::domainIntegrate(qDecay_).value()*dimensionedScalar("", dimPower, 1);

    Q_ = Qprompt_ + Qdecay_;
}
