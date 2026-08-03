{
    forAll(flux_, i)
    {
        // Temp zero fields
        volScalarField tSNS = zeroField_*dimensionedScalar("", dimless/dimVolume/dimTime, 0);
        volScalarField tFNS = zeroField_*dimensionedScalar("", dimless/dimVolume/dimTime, 0);
        volScalarField tDNS = zeroField_*dimensionedScalar("", dimless/dimVolume/dimTime, 0);

        // Evaluate neutron sources based on current fluxes and cross sections
        forAll(flux_, j)
        {
            if (j != i)
            {
                tSNS += Ss_[nEnergyGroups_*j + i]*flux_[j];
                tFNS += chiPrompt_[i]*(1 - precBetaTot_)*Nu_[j]*Sf_[j]*flux_[j];
            }
        }

        forAll(prec_, j)
        {
            tDNS += chiDelayed_[i]*precLambda_[j]*prec_[j];
        }

        // Repopulate pointer lists
        neutronSource_[i] = tSNS + tFNS + tDNS;
    }
}
