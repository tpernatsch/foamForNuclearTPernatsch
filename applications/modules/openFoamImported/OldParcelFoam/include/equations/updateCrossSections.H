{
    // Evaluate density and temperature fields
    rhok_ = 1.0 - beta_*(T_ - TRef_);
    rhoNorm_ = rhok_/(1.0 - beta_*(TRefXS_ - TRef_));
    logT_ = log(max(0.5, T_/TRefXS_));

    // Evaluate temperature- and density-dependent nuclear properties
    forAll(flux_, i)
    {
        D_[i] = (D_ref_[i] + DopplerControl_*alphaD_ref_[i]*logT_)/rhoNorm_;
        Sa_[i] = (Sa_ref_[i] + DopplerControl_*alphaSa_ref_[i]*logT_)*rhoNorm_;
        Sf_[i] = (1/Keff_)*(Sf_ref_[i] + DopplerControl_*alphaSf_ref_[i]*logT_)*rhoNorm_;
    }

    forAll(Ss_, i)
    {
        Ss_[i] = (Ss_ref_[i] + DopplerControl_*alphaSs_ref_[i]*logT_)*rhoNorm_;
    }

    // Evaluate net neutron production operator
    forAll(flux_, i)
    {
        // Total removal cross-section of i-th group due to scattering to other groups
        volScalarField Ssi = 0*Ss_[0];

        forAll(flux_, j)
        {
            if (j != i)
            {
                Ssi += Ss_[nEnergyGroups_*i + j];
            }
        }

        // Neutrons produced from fissions in i-th group minus neutrons removed
        M_[i] = chiPrompt_[i]*(1 - precBetaTot_)*Nu_[i]*Sf_[i] - (Sa_[i] + Ssi);
    }
}
