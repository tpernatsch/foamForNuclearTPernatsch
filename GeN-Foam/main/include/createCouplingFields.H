//- Determine if thermalHydraulics is onePhase or twoPhase by counting fluid
//  objects
HashTable<fluid*> fluids(fluidMesh.lookupClass<fluid>());
int nFluids(fluids.size());

//- Non-const pointers to non-const fluid objects. The pointers are non-const
//  as they need to be set here. If only one fluid is present, no big deal. If
//  two are present, determine which is the vapour and which is the liquid.
//  This is based by comparing enthalpies of formation as specified in the
//  thermophysicalProperties dictionaries. The liquidFuel option (e.g. for
//  modelling MSRs) is only meaningful for single-phase simulations. Ever heard
//  of boiling MSRs? Yeah, me neither...
fluid* liquid(nullptr);
fluid* vapour(nullptr);
if (nFluids == 1)
{
    liquid = fluids[fluids.toc()[0]];

    //- If the liquidFuel option is on, the powerDensity will be inserted in
    //  the liquid energy equation automatically by the thermalHydraulic
    //  class when doing correctEnergy. However, the actual fluid powerDensity
    //  field needs to be created (i.e., the autoPtr to said field within the
    //  fluid class). Why did I use an autoPtr, and not just a regular
    //  volScalarField? Well, that would force the fluid class to always create
    //  a power density field (also for the vapour), and that is a useless
    //  waste of RAM and disk memory (the field needs to be written)
    if (multiphysics.liquidFuel())
    {
        liquid->powerDensityPtr().reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName("powerDensity", liquid->name()),
                    fluidMesh.time().timeName(),
                    fluidMesh.time(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                fluidMesh,
                dimensionedScalar("", dimPower/dimVol, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }
}
else
{   
    if (multiphysics.liquidFuel())
    {
        FatalErrorInFunction
        << "Simulation of twoPhase scenarios involving liquidFuel is "
        << "unsupported"
        << exit(FatalError);
    }
    
    fluid* const fluid0(fluids[fluids.toc()[0]]);
    fluid* const fluid1(fluids[fluids.toc()[1]]);

    //- Use average enthalpies of formation for comparison. Doesn't really
    //  matter unless one has a space-dependent enthalpy of formation, which
    //  I have not seen so far in OpenFOAM
    scalar hc0(gAverage(fluid0->thermo().hc()()));
    scalar hc1(gAverage(fluid1->thermo().hc()()));

    if (hc1 > hc0)
    {
        liquid = fluid0;
        vapour = fluid1;
    }
    else if (hc0 > hc1)
    {
        vapour = fluid0;
        liquid = fluid1;
    }
    else
    {
        FatalErrorInFunction
            << "When running in twoPhase mode, vapour and liquid are "
            << "distinguished by their enthalpy of formation, i.e. the Hf "
            << "keyword specified in the thermophysicalProperties of each "
            << "fluid. " << endl
            << "The fluid with the higher Hf will be the vapour, and the "
            << "latent heat of evaporation will be difference between the two."
            << "Currently, both fluids appear to have the same Hf. "
            << exit(FatalError);
    }
}

//- Get ref to passive structure temperature
const volScalarField& TStruct
(
    (*fluidMesh.findObject<volScalarField>("T.passiveStructure"))
);

//- Create pointers to fuel and cladding average temperatures. These might not
//  exist, it depends on the existence of the nuclearFuelPin models. Hence, it
//  is more convenient to have these as ptrs rather than references
volScalarField* const TFuelAv
(
    (fluidMesh.foundObject<volScalarField>("Tfav.nuclearFuelPin"))  ?
    fluidMesh.findObject<volScalarField>("Tfav.nuclearFuelPin") :
    new volScalarField
    (
        IOobject
        (
            "nullTFuelAv0",
            fluidMesh.time().timeName(),
            fluidMesh.time()
        ),
        fluidMesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
);
volScalarField* const TCladAv
(
    (fluidMesh.foundObject<volScalarField>("Tcav.nuclearFuelPin"))  ?
    fluidMesh.findObject<volScalarField>("Tcav.nuclearFuelPin") :
    new volScalarField
    (
        IOobject
        (
            "nullTCladAv0",
            fluidMesh.time().timeName(),
            fluidMesh.time()
        ),
        fluidMesh,
        dimensionedScalar("", dimTemperature, 0.0),
        zeroGradientFvPatchScalarField::typeName
    )
);

//- Set the actual volFuelPower field reference, that gets modfied by
//  the neutronics for coupling purposes. If doing liquidFuel, the volFuelPower
//  field is the liquid powerDensity. Otherwise, if a nuclearFuelPin-related
//  powerDensity field exists (which is true if at least one nuclearFuelPin
//  powerModel figures in phaseProperties), use that one. Otherwise again, if
//  no powerDensity field exists, create a null one. The volFuelPower field is 
//  a ptr and not a reference as references cannot be initialized with rvalues,
//  such as is the new volScalarField (i.e., I would to create it outside the 
//  ptr initialization, and that again would happen always and just waste 
//  memory if the other powerDensity fields exist). I do realize these are 
//  laughably small amounts, but hey, whatever helps
volScalarField* const volFuelPower
(
    (multiphysics.liquidFuel()) ? &(liquid->powerDensity()) :
    (
        (fluidMesh.foundObject<volScalarField>("powerDensity.nuclearFuelPin")) 
        ?
        fluidMesh.findObject<volScalarField>("powerDensity.nuclearFuelPin") :
        new volScalarField
        (
            IOobject
            (
                "nullPowerDensity",
                fluidMesh.time().timeName(),
                fluidMesh.time()
            ),
            fluidMesh,
            dimensionedScalar("", dimless, 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    )
);

//- Silence "unused variable" warnings, I actually do use them and have no idea
//  why the compiler does not "see" them
(void) liquid;
(void) vapour;
(void) TStruct;
(void) TFuelAv;
(void) TCladAv;
(void) volFuelPower;

//- Then, create the actual fields for the coupling. These are created over the
//  fluidMeshNB mesh, which is a baffle-less mesh that is created if the 
//  removeBaffles option is set in controlDict. This is generally required when
//  doing multi-physics on meshes that contain baffles, as baffles are known to
//  mess up mesh-to-mesh projections. If the removeBaffles options is not
//  specified, fluidMeshNB defaults to fluidMesh (see createFluidMesh.H)
volScalarField bafflelessTCool
(
    IOobject
    (
        "bafflelessTCoolant",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessTFuelAv
(
    IOobject
    (
        "bafflelessTFuelAv",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("",  dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessTCladAv
(
    IOobject
    (
        "bafflelessTavClad",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("",  dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessTStruct
(
    IOobject
    (
        "bafflelessTStruct",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimTemperature, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessAlpha
(
    IOobject
    (
        "bafflelessAlpha",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimless, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField bafflelessRho
(
    IOobject
    (
        "bafflelessRho",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimDensity, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField bafflelessU
(
    IOobject
    (
        "bafflelessU",
        runTime.timeName(),
        fluidMesh
    ),
    fluidMeshNB,
    dimensionedVector("", dimVelocity, vector::zero),
    zeroGradientFvPatchScalarField::typeName
);

//- This last one is used in place of mu, alphat, as I noticed they are not
//  used at all by any neutronics subsolver...
volScalarField bafflelessNullField
(
    IOobject
    (
        "bafflelessNullField",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimless, 0.0)
);
