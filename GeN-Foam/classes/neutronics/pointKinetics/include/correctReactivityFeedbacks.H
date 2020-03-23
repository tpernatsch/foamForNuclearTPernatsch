//- Create field for weighting. This is currently done on the oneGroupFlux,
//  rather than on the energy-dependent adjoint. It will be changed in the
//  future
volScalarField oneGroupFluxSqr(sqr(oneGroupFlux_));

scalar domainIntegratedFluxSqrFuel
(
    fvc::domainIntegrate(oneGroupFluxSqr*fuelFeedbackCellField_).value()
);
scalar domainIntegratedFluxSqrCoolant
(
    fvc::domainIntegrate(oneGroupFluxSqr*coolantFeedbackCellField_).value()
);

//- Compute average perturbed values
scalar TFuel
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TFuel_*fuelFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrFuel
);
scalar TCladding
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TCladding_*fuelFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrFuel
);
scalar TCoolant
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TCoolant_*coolantFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrCoolant
);
scalar rhoCoolant
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*rhoCoolant_*coolantFeedbackCellField_
    ).value()/
    domainIntegratedFluxSqrCoolant
);
scalar TStructures
(
    fvc::domainIntegrate
    (
        oneGroupFluxSqr*TStructures_*structuresFeedbackCellField_
    ).value()/
    fvc::domainIntegrate(oneGroupFluxSqr*structuresFeedbackCellField_).value()
);

//- Not weighted by flux
scalar TDriveline
    (
        fvc::domainIntegrate
        (
            TStructures_*drivelineFeedbackCellField_
        ).value()/
        max
        (
            fvc::domainIntegrate(drivelineFeedbackCellField_).value(),
            SMALL
        )
    );

scalar drivelineReactivityValue(0);    //- Name conflict with the 
                                        //  drivelineReactivity function
scalar GEMReactivityValue(0); //- Name conflict with the GEMReactivity function
scalar DopplerReactivity(0);
scalar TFuelReactivity(0);
scalar TCladdingReactivity(0);
scalar TCoolantReactivity(0);
scalar rhoCoolantReactivity(0);
scalar TStructuresReactivity(0);

//- If the user is starting from a steady state, I need to compute
//  the ref values. But the ref values require knowledge of the 
//  fields in the fluidMesh, and the neutronics only gets them
//  AFTER construction, it cannot ask for them before (this should be
//  changed in the rest of the GeN-Foam, outside this class, somehow).
//  The best I can do for know is to assume that the reactor is at steady 
//  state at time-step number 1, rather than the starting time-step.
if  
(
        mesh_.time().timeIndex()
    -   mesh_.time().startTimeIndex() == 1 
)
{
    TFuelRef_ = resetTFuelRef_ ? TFuel : TFuelRef_;
    TCladdingRef_ = resetTCladdingRef_ ? TCladding : TCladdingRef_;
    TCoolantRef_ = resetTCoolantRef_ ? TCoolant : TCoolantRef_;
    rhoCoolantRef_ = resetRhoCoolantRef_ ? rhoCoolant : rhoCoolantRef_;
    TStructuresRef_ = resetTStructuresRef_ ? TStructures : TStructuresRef_;
    TDrivelineRef_ = TDriveline; //- I'm rushing through this

    //- Set reference flow value for GEM modelling
    flowRef_ *= 0;
    forAll(flowFaces_, i)
    {
        flowRef_ += phi_[flowFaces_[i]];
    }
    //- Sync integral flow value across processors
    reduce(flowRef_, sumOp<scalar>());
    flowRef_ = std::max(std::abs(flowRef_), 1e-69);
}
else
{
    drivelineReactivityValue = drivelineReactivity(); 
    GEMReactivityValue = GEMReactivity();
    DopplerReactivity = -coeffFastDoppler_*Foam::log(TFuel/TFuelRef_);
    TFuelReactivity = coeffTFuel_*(TFuel - TFuelRef_);
    TCladdingReactivity = coeffTCladding_*(TCladding - TCladdingRef_);
    TCoolantReactivity = coeffTCoolant_*(TCoolant - TCoolantRef_);
    rhoCoolantReactivity = coeffRhoCoolant_*(rhoCoolant - rhoCoolantRef_);
    TStructuresReactivity = coeffTStructures_*(TStructures - TStructuresRef_);

    totalReactivity_ = 
            externalReactivity_ 
        +   DopplerReactivity
        +   TFuelReactivity
        +   TCladdingReactivity
        +   TCoolantReactivity
        +   rhoCoolantReactivity
        +   TStructuresReactivity
        +   GEMReactivityValue
        +   drivelineReactivityValue;

    reactivityValues_["total"] = totalReactivity_;
    reactivityValues_["external"] = externalReactivity_;
    reactivityValues_["Doppler"] = DopplerReactivity;
    reactivityValues_["TFuel"] = TFuelReactivity;
    reactivityValues_["TCladding"] = TCladdingReactivity;
    reactivityValues_["TCoolant"] = TCoolantReactivity;
    reactivityValues_["rhoCoolant"] = rhoCoolantReactivity;
    reactivityValues_["TStructures"] = TStructuresReactivity;
    reactivityValues_["GEM"] = GEMReactivityValue;
    reactivityValues_["driveline"] = drivelineReactivityValue;
}
