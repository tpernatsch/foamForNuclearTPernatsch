//- Temperature and density feedbacks

scalar DopplerReactivity = 0.0;
if (fastNeutrons_)
{
    DopplerReactivity = -coeffFastDoppler_*
        Foam::log(max(TFuelValue, SMALL) / max(TFuelRef_, SMALL));
} 
else // For thermal spectrum
{
    DopplerReactivity = coeffFastDoppler_* (sqrt(TFuelValue) - sqrt(TFuelRef_));
}

scalar TFuelReactivity =
    coeffTFuel_*(TFuelValue - TFuelRef_);

scalar TCladReactivity =
    coeffTClad_*(TCladValue - TCladRef_);

scalar TCoolReactivity =
    coeffTCool_*(TCoolValue - TCoolRef_);

scalar rhoCoolReactivity =
    coeffRhoCool_*(rhoCoolValue - rhoCoolRef_);

scalar TStructReactivity =
    coeffTStruct_*(TStructValue - TStructRef_);

scalar TStructMechReactivity =
    coeffTStructMech_*(TStructMechValue - TStructMechRef_);

//- Driveline reactivity
scalar drivelineExpValue
(
    coeffDrivelineExp_*(TDrivelineValue - TDrivelineRef_)
);
scalar drivelineReactivity(calcDrivelineReactivity(drivelineExpValue));

//- GEM reactivity
Pair<scalar> GEMLevelAndReactivity(calcGEMLevelAndReactivity());
scalar GEMSodiumLevel(GEMLevelAndReactivity.first());
scalar GEMReactivity(GEMLevelAndReactivity.second());

(void) GEMSodiumLevel;

//- External reactivity
scalar t(mesh_.time().timeOutputValue());
scalar externalReactivity = externalReactivityTimeProfile_.value(t);

//- FMU
scalar externalReactivityFromFMU(0.0);
#ifdef isCommDataLayerIncluded

word externalReactivityKeyFromFMU("externalReactivityNameFromFMU");
if (nuclearData_.found(externalReactivityKeyFromFMU))
{
    const word externalReactivityNameFromFMU = nuclearData_.get<word>(externalReactivityKeyFromFMU);

    // Communicating with the FMU
    const Time& runTime = this->db().time();
    commDataLayer& data = commDataLayer::New(runTime);
    externalReactivityFromFMU = data.getObj<scalar>(
        externalReactivityNameFromFMU,
        commDataLayer::causality::in
    );
}

#endif // isCommDataLayerIncluded


//- Boron reactivity
scalar boronReactivity = boronReactivityTimeProfile_.value(t);

//- Update total reactivity
totalReactivity_ =
        precEquilibriumReactivity_
    +   externalReactivity
    +   externalReactivityFromFMU
    +   DopplerReactivity
    +   TFuelReactivity
    +   TCladReactivity
    +   TCoolReactivity
    +   rhoCoolReactivity
    +   TStructReactivity
    +   TStructMechReactivity
    +   drivelineReactivity
    +   GEMReactivity
    +   boronReactivity;
