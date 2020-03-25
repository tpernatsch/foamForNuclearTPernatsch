scalar DopplerReactivity = 
    -coeffFastDoppler_*Foam::log(TFuelValue/max(TFuelRef_, SMALL));

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

totalReactivity_ = 
        externalReactivity_ 
    +   DopplerReactivity
    +   TFuelReactivity
    +   TCladReactivity
    +   TCoolReactivity
    +   rhoCoolReactivity
    +   TStructReactivity;