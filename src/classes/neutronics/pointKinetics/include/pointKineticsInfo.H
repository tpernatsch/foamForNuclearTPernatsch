Info<< "    totalPower = " << power_ << " W" << nl
    << "    -> fission = " << fissionPower_ << " W" << nl
    << "    -> decay   = " << decayPower_ << " W" << nl
    << "    totalReactivity   = " << (totalReactivity_*1e5)<< " pcm" << nl
    << "    -> extReactivity  = " << (externalReactivity*1e5)<< " pcm" << nl
    << "    -> Doppler (fast) = " << (DopplerReactivity*1e5)<< " pcm" << nl
    << "    -> TFuel          = " << (TFuelReactivity*1e5) << " pcm" << nl
    << "    -> TClad          = " << (TCladReactivity*1e5) << " pcm" << nl
    << "    -> TCool          = " << (TCoolReactivity*1e5) << " pcm" << nl
    << "    -> rhoCool        = " << (rhoCoolReactivity*1e5) << " pcm" << nl
    << "    -> TStruct        = " << (TStructReactivity*1e5) << " pcm" << nl
    << "    -> TStructMech    = " << (TStructMechReactivity*1e5) << " pcm" << nl
    << "    -> driveline      = " << (drivelineReactivity*1e5) << " pcm"<< nl
    << "    -> Boron          = " << (boronReactivity*1e5) << " pcm"
    << endl;

if (externalSourceNeutronics_)
{
    Info<< "    -> subcriticalIdx = " << (subcriticalIndex_*1e5) << " pcm" << nl
        << "    extSrcPower         = " << externalSourcePower_ << " W" << nl
        << "    extSrcBeamPower     = " << externalSourceBeamPower_ << " W" << nl
        << "    extSrcBeamIntensity = " << externalSourceBeamIntensity_*1.602176634e-19 << " A" << nl
        << "    extSrcModulation    = " << externalSourceModulation_
        << endl;
}

if (liquidFuel_)
{
    Info<< "    Initial reactivity loss due to delayed neutron precursors movement = "
        << (precEquilibriumReactivity_ *1e5) << " pcm"
        << nl
        << "    Estimated Beta including DNP recirculation = "
        << (liquidFuelBeta_ *1e5) << " pcm"
        << endl;
}

if (this->get<bool>("GEM"))
{
	Info<< "    -> GEM            = " << (GEMReactivity*1e5) << " pcm" << endl;
}

Info<< "    TFuel = " << TFuelValue << " K" << nl
    << "    TClad = " << TCladValue << " K" << nl
    << "    TCool = " << TCoolValue << " K" << nl
    << "    rhoCool = " << rhoCoolValue << " kg/m3" << nl
    << "    TStruct = " << TStructValue << " K" << nl
    << "    TStructMech = " << TStructMechValue << " K" << nl
    << "    driveline (T, exp) = " << TDrivelineValue << " K, " << drivelineExpValue << " m"
    << endl;

if (this->get<bool>("GEM"))
	Info<< "    GEMSodiumLevel = " << GEMSodiumLevel << " m" << endl;
