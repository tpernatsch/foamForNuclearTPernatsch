//- Brief note on use of the const keyword with pointers:
/*
    Type* ptr;  // ptr CAN be reassigned and points to an object of type Type
                // that CAN be modified via ptr

    const Type* ptr; // ptr CAN be reassigned and points to an object of type
                // Type that CANNOT be modified via ptr

    Type* const ptr; // ptr CANNOT be reassigned and points to an object of 
                // type Type that CAN be modified via ptr

    const Type* const ptr // ptr CANNOT be reassigned and points to an object 
                // of type Type that CANNOT be modified via ptr
*/

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
                    liquid->db(),
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

//- The powerDensity field is created exclusively by the nuclearFuelPin,
//  so if it does exists, so does the nuclearFuelPin model
const bool nuclearFuelPin
(
    (fluidMesh.foundObject<volScalarField>("powerDensity.nuclearFuelPin")) 
);

//- Create pointers to fuel and cladding average temperatures. These might not
//  exist, it depends on the existence of the nuclearFuelPin models. Hence, it
//  is more convenient to have these as ptrs rather than references
volScalarField* const TFuelAv
(
    (nuclearFuelPin)  ?
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
    (nuclearFuelPin)  ?
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
//  such as is the new volScalarField (i.e., I would nned to create it outside  
//  the ptr initialization, and that again would happen always and just waste 
//  memory if the other powerDensity fields exist). I do realize these are 
//  laughably small of memory to bitch about, but hey, whatever helps
volScalarField* const volFuelPower
(
    (multiphysics.liquidFuel()) ? &(liquid->powerDensity()) :
    (
        (nuclearFuelPin) 
        ?
        fluidMesh.findObject<volScalarField>("powerDensity.nuclearFuelPin") :
        new volScalarField
        (
            IOobject
            (
                "powerDensity.nuclearFuelPin",
                fluidMesh.time().timeName(),
                liquid->db(),
                IOobject::READ_IF_PRESENT,
                (
                    (
                        multiphysics.solveEnergy()
                    or  multiphysics.solveNeutronics()
                    ) ? IOobject::AUTO_WRITE : IOobject::NO_WRITE
                )
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
        "bafflelessTCool",
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
        "bafflelessTCladAv",
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

volScalarField bafflelessRhoCool
(
    IOobject
    (
        "bafflelessRhoCool",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimDensity, SMALL),
    zeroGradientFvPatchScalarField::typeName
);

//- Turbulent thermal diffusivity of liquid. For some reason, thermal
//  diffusivity as returned by thermo() is always multiplied by the density and
//  is not in m2/s, rather, kg/m/s (== Pa*s). It's a standard OpenFOAM thing,
//  not a bug, yet still...
volScalarField bafflelessAlphat
(
    IOobject
    (
        "bafflelessAlphat",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimArea/dimTime*dimDensity, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

//- Molecular viscosity of liquid
volScalarField bafflelessMu
(
    IOobject
    (
        "bafflelessMu",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedScalar("", dimPressure*dimTime, 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volVectorField bafflelessU
(
    IOobject
    (
        "bafflelessU",
        runTime.timeName(),
        fluidMeshNB
    ),
    fluidMeshNB,
    dimensionedVector("", dimVelocity, vector::zero),
    zeroGradientFvPatchScalarField::typeName
);

//- This field is sort of an exception required when using pointKinetics.
//  PointKinetics by itself simply rescales all neutronic fields according to
//  the predicted power change. However, if no volFuelPower exists on the
//  neutroMesh, then look for volFuelPower on the fluidMesh. If it exists, copy
//  it to a bafflelessVolFuelPower that will be interpolated to the
//  volFuelPower_ on the neutroMesh of pointKinetics in the 
//  getCouplingFieldRefs function. Since pointKinetics is also just about
//  rescaling an initial field, this bafflelessVolFuelPower will not need to be
//  set every time, only the first. Nedless to say, pointKinetics can function
//  perfectly well with a 0 volFuelPower, it simply won't be re-scaling 
//  anything (but integral power values will be predicted just fine)
volScalarField* bafflelessVolFuelPower = nullptr;

//- This is only required when point kinetics model Gas Expansion Modules (i.e.
//  a safety feature specific to one SFR, the FFTF. I hate this level of hard-
//  coding but I don't have alternatives for now). It is not baffleless as this
//  field is not interpolated, it is required on the original fluidMesh by the
//  pointKinetics model, yet it needs to be inserted in the fluidMeshNB 
//  objectRegistry, or the pointKinetics model won't be able to look it up.
//  phiFMNBReg is the phi but checked in the the fluidMeshNB registry, while
//  phiFMReg is the phi in the fluidMesh registry (i.e., the original one).
//  both are however defined on the fluidMesh mesh anyway
surfaceScalarField* phiFMNBReg = nullptr;
const surfaceScalarField& phiFMReg
(
    (fluidMesh.lookupObject<surfaceScalarField>("phi"))
);

if (neutronics.type() == "pointKinetics")
{
    bafflelessVolFuelPower = 
        new volScalarField
        (
            IOobject
            (
                "bafflelessVolFuelPower",
                runTime.timeName(),
                fluidMeshNB
            ),
            fluidMeshNB,
            dimensionedScalar("", dimPower/dimDensity, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );
    bafflelessVolFuelPower->primitiveFieldRef() = 
        volFuelPower->primitiveField();
    bafflelessVolFuelPower->correctBoundaryConditions();
    
    //- Phi only required by GEM modelling
    if (neutronics.get<bool>("GEM"))
    {
        phiFMNBReg =
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidMeshNB
                ),
                fluidMesh,
                dimensionedScalar("", dimVolume/dimTime, 0.0)
            );
        *phiFMNBReg = phiFMReg;
    }
}

//- Update the above fields so that, if any neutronic classes rely on these
//  fields for further local field initialization (e.g. pointKinetics), correct
//  initial conditions are set
if (multiphysics.solveNeutronics() or multiphysics.solveThermoMechanics())
{
    #include "correctCouplingFields.H"
}

//- Set local neutronics and thermoMechanics references to fields they require
//  for coupling. These are all obtained via lookups. Additionally, the
//  neutronics->thermalHydraulic interpolator is passed as some classes might
//  require to project some fields right after construction and before the 
//  start of the simulation e.g. pointKinetics
if (multiphysics.solveNeutronics())
{
    neutronics.getCouplingFieldRefs
    (
        fluidMeshNB,
        neutroToFluid
    );
}
if (multiphysics.solveThermoMechanics())
{
    thermoMechanics.getCouplingFieldRefs
    (
        fluidMeshNB,
        mechToFluid
    );
}
