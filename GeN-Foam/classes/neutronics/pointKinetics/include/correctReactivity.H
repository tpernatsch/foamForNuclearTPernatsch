//- Temperature and density feedbacks
scalar DopplerReactivity =
    -coeffFastDoppler_*
    Foam::log(max(TFuelValue, SMALL)/max(TFuelRef_, SMALL));

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

//- Boron reactivity
scalar boronReactivity = boronReactivityTimeProfile_.value(t);

//- Update total reactivity
totalReactivity_ =
        precEquilibriumReactivity_
    +   externalReactivity
    +   DopplerReactivity
    +   TFuelReactivity
    +   TCladReactivity
    +   TCoolReactivity
    +   rhoCoolReactivity
    +   TStructReactivity
    +   drivelineReactivity
    +   GEMReactivity
    +   boronReactivity;
