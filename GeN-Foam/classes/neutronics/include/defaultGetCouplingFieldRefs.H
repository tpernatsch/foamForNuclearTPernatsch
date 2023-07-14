//- Field names must reflect those defined in createCouplingFields.H
TFuelOrig_ = 
    srcTH.findObject<volScalarField>("bafflelessTFuelAv");
TCladOrig_ = 
    srcTH.findObject<volScalarField>("bafflelessTCladAv");
TCoolOrig_ = 
    srcTH.findObject<volScalarField>("bafflelessTCool");
rhoCoolOrig_ = 
    srcTH.findObject<volScalarField>("bafflelessRhoCool");
if (liquidFuel_)
{
    UOrig_ = 
        srcTH.findObject<volVectorField>("bafflelessU");
    alphaOrig_ = 
        srcTH.findObject<volScalarField>("bafflelessAlpha");
    alphatOrig_ = 
        srcTH.findObject<volScalarField>("bafflelessAlphat");
    muOrig_ =
        srcTH.findObject<volScalarField>("bafflelessMu");
}
else
{
    UOrig_ = nullptr;
    alphaOrig_ = nullptr;
    alphatOrig_ = nullptr;
    muOrig_ = nullptr;
}

//- Get T from TM solver
TStructMechOrig_ = 
    srcTM.findObject<volScalarField>("TStruct");

//- Initialize mapped fields
this->interpolateCouplingFields(neutroToFluid,neutroToMech);

//- Initialize XS correctly
xs_.correct(TFuel_, TClad_, rhoCool_, TCool_, disp_, TStructMech_);
