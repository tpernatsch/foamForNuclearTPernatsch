// This include modifies the dictionaries so that GeN-Foam V2 becomes the equivalent of a GeN-Foam V1 case


const bool solveThermalHydraulics(runTime.controlDict().getOrDefault<bool>("solveFluidMechanics", false));
const bool solveEnergy(runTime.controlDict().getOrDefault<bool>("solveEnergy", false));
const bool solveNeutronics(runTime.controlDict().getOrDefault<bool>("solveNeutronics", false));
const bool solveThermoMechanics(runTime.controlDict().getOrDefault<bool>("solveThermalMechanics", false));
const bool liquidFuel(runTime.controlDict().getOrDefault<bool>("liquidFuel", false));
bool PK(false);
bool linkedFuel(false);

// First: check whether simulation is tightly coupled or not, to see whether a multiPhysicsSolver has to be created.


IOdictionary multiRegionCouplingDict
(
    IOobject
    (
        "multiRegionCouplingDict",
        runTime.constant(),
        runTime.db(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    )
);



IOdictionary fvSolutionDict
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        runTime.db(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

bool tightlyCoupled(fvSolutionDict.get<bool>("tightlyCoupled"));

dictionary multiPhysicsSolver;
dictionary subSolvers;

IOdictionary controlDict
(
    IOobject
    (
        "controlDict",
        runTime.system(),
        runTime.db(),
        IOobject::MUST_READ,
        IOobject ::NO_WRITE
    )
);

dictionary regionSolversDict;

if(controlDict.found("regionSolvers"))
{
    controlDict.remove("regionSolvers");
}

controlDict.regIOobject::write();

if(solveThermalHydraulics or solveEnergy)
{
    const IOdictionary phaseProperties
    (
        IOobject
        (
            "phaseProperties",
            "constant/fluidRegion",
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word fluidSolver(phaseProperties.get<word>("thermalHydraulicsType"));

    if (!tightlyCoupled)
        regionSolversDict.add("fluidRegion", fluidSolver);
    else
        subSolvers.add("fluidRegion", fluidSolver);


    if (solveEnergy)
    {
        IOdictionary fvSolTH
        (
            IOobject
            (
                "fvSolution",
                runTime.system(),
                "fluidRegion",
                runTime.db(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            )
        );

        fvSolTH.subDict("PIMPLE").add("solveEnergy", "true");
        if (!solveThermalHydraulics)
        {
            fvSolTH.subDict("PIMPLE").add("solveFluidMechanics", "false");
        }
        fvSolTH.regIOobject::write();
    }


}

if (solveEnergy)
{
    IOdictionary fvSolTH
    (
        IOobject
        (
            "fvSolution",
            runTime.system(),
            "fluidRegion",
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        )
    );

    fvSolTH.subDict("PIMPLE").add("solveEnergy", "true");
    if (!solveThermalHydraulics)
    {
        fvSolTH.subDict("PIMPLE").add("solveFluidMechanics", "false");
    }
    fvSolTH.regIOobject::write();

}

if (solveNeutronics)
{
    const IOdictionary neutronicsProperties
    (
        IOobject
        (
            "neutronicsProperties",
            "constant/neutroRegion",
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word neutroSolver(neutronicsProperties.get<word>("model"));

    if (neutroSolver == "pointKinetics")
        PK = true;

    if (!tightlyCoupled)
        regionSolversDict.add("neutroRegion", neutroSolver);
    else
        subSolvers.add("neutroRegion", neutroSolver);
}

if (solveThermoMechanics)
{
    const IOdictionary thermoMechProperties
    (
        IOobject
        (
            "thermoMechanicalProperties",
            "constant/thermoMechanicalRegion",
            runTime.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const word mechSolver(thermoMechProperties.get<word>("model"));
    linkedFuel = thermoMechProperties.getOrDefault<bool>("linkedFuel", false);

    if (!tightlyCoupled)
        regionSolversDict.add("thermoMechanicalRegion", mechSolver);
    else
        subSolvers.add("thermoMechanicalRegion", mechSolver);
}


if (tightlyCoupled)
    multiPhysicsSolver.add("solvers",subSolvers);




//- Now add the removeBafflesOption

bool removeFluidBaffles(controlDict.getOrDefault<bool>("removeFluidBaffles", false));
dictionary removeBafflesDict;
controlDict.remove("removeFluidBaffles");
controlDict.remove("removeBaffles");
removeBafflesDict.add("fluidRegion", removeFluidBaffles);
removeBafflesDict.add("neutroRegion", false);
removeBafflesDict.add("thermoMechanicalRegion", false);

Info<< removeBafflesDict << endl;
controlDict.add("removeBaffles", removeBafflesDict);



// - Now create and write the multiRegionCouplingDict
if (tightlyCoupled)
{
    multiPhysicsSolver.add("maxIterations",  fvSolutionDict.get<scalar>("maxTimeStepIterations"));
    multiPhysicsSolver.add("minResidual",  fvSolutionDict.get<scalar>("timeStepResidual"));
    dictionary solverName;
    solverName.add("tightLoop", multiPhysicsSolver);
    multiRegionCouplingDict.add("multiPhysicsSolvers", solverName);
    regionSolversDict.add("tightLoop", "multiPhysicsSolver");
}

controlDict.add("regionSolvers", regionSolversDict);
// Add outer correctors of the TH cycle

IOdictionary fvSolTH
(
    IOobject
    (
        "fvSolution",
        runTime.system(),
        "fluidRegion",
        runTime.db(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    )
);

fvSolTH.subDict("PIMPLE")
    .add("nOuterCorrectors", fvSolutionDict.get<scalar>("nOuterCorrectors"));

fvSolTH.regIOobject::write();

// - Add the mappings list
dictionary mappings;

if (solveThermalHydraulics or solveEnergy)
{
    dictionary THMappingDict;

    if (solveNeutronics)
    {
        dictionary neuToTH;
        List<word> tgtNames(2);
        tgtNames[0] = "powerDensityNeutronics";
        tgtNames[1] = "powerDensityNeutronicsToLiquid";
        List<word> fieldTypes(2);
        fieldTypes[0] = "scalar";
        fieldTypes[1] = "scalar";
        List<word> fieldNames(2);

        if (!liquidFuel)
        {
            fieldNames[0] = "powerDensity";
            fieldNames[1] = "secondaryPowerDensity";
        }
        else
        {
            fieldNames[1] = "powerDensity";
            fieldNames[0] = "secondaryPowerDensity";
        }
        neuToTH.add("sourceFields", fieldNames);
        neuToTH.add("targetFields", tgtNames);
        neuToTH.add("fieldTypes", fieldTypes);

        THMappingDict.add("neutroRegion", neuToTH);
    }

    mappings.add("fluidRegion", THMappingDict);
}


if (solveNeutronics)
{
    dictionary NeuMappingDict;

    if (solveEnergy or solveThermalHydraulics)
    {
        dictionary THToNeu;
        if(!liquidFuel)
        {
            List<word> tgtNames(5);
            tgtNames[0] = "TCool";
            tgtNames[1] = "rhoCool";
            tgtNames[2] = "TFuel";
            tgtNames[3] = "TClad";
            tgtNames[4] = "TStructMech";
            List<word> fieldTypes(5);
            fieldTypes[0] = "scalar";
            fieldTypes[1] = "scalar";
            fieldTypes[2] = "scalar";
            fieldTypes[3] = "scalar";
            fieldTypes[4] = "scalar";
            List<word> fieldNames(5);
            fieldNames[0] = "T";
            fieldNames[1] = "thermo:rho";
            fieldNames[2] = "T.fuelAvForNeutronics";
            fieldNames[3] = "T.cladAvForNeutronics";
            fieldNames[4] = "T.passiveStructure";

            THToNeu.add("sourceFields", fieldNames);
            THToNeu.add("targetFields", tgtNames);
            THToNeu.add("fieldTypes", fieldTypes);


            NeuMappingDict.add("fluidRegion", THToNeu);
        }

        else
        {
            List<word> tgtNames(9);
            tgtNames[0] = "TCool";
            tgtNames[1] = "rhoCool";
            tgtNames[2] = "TFuel";
            tgtNames[3] = "TClad";
            tgtNames[4] = "TStructMech";
            tgtNames[5] = "U";
            tgtNames[6] = "mu";
            tgtNames[7] = "alpha";
            tgtNames[8] = "alphat";
            List<word> fieldTypes(9);
            fieldTypes[0]= "scalar";
            fieldTypes[1]= "scalar";
            fieldTypes[2]= "scalar";
            fieldTypes[3]= "scalar";
            fieldTypes[4]= "scalar";
            fieldTypes[5]= "vector";
            fieldTypes[6]= "scalar";
            fieldTypes[7]= "scalar";
            fieldTypes[8]= "scalar";
            List<word> fieldNames(9);
            fieldNames[0]= "T";
            fieldNames[1]= "thermo:rho";
            fieldNames[2]= "T";
            fieldNames[3]= "T.cladAvForNeutronics";
            fieldNames[4] = "T.passiveStructure";
            fieldNames[5] = "U";
            fieldNames[6] = "mu";
            fieldNames[7] = "alpha";
            fieldNames[8] = "thermo:alpha";

            THToNeu.add("sourceFields", fieldNames);
            THToNeu.add("targetFields", tgtNames);
            THToNeu.add("fieldTypes", fieldTypes);
        }

        NeuMappingDict.add("fluidRegion", THToNeu);

    }

    if (solveThermoMechanics and PK)
    {
        dictionary TMToNeu;

        List<word> tgtNames(1);
        tgtNames[0] = "TStructMech";
        List<word> fieldTypes(1);
        fieldTypes[0] = "scalar";
        List<word> fieldNames(1);
        fieldNames[0] = "TStruct";

        TMToNeu.add("sourceFields", fieldNames);
        TMToNeu.add("targetFields", tgtNames);
        TMToNeu.add("fieldTypes", fieldTypes);

        NeuMappingDict.add("thermoMechanicalRegion", TMToNeu);
    }

    if (solveThermoMechanics and !PK)
    {
        dictionary TMToNeu;
        List<word> tgtNames(1);
        tgtNames[0] = "disp";
        List<word> fieldTypes(1);
        fieldTypes[0] = "vector";
        List<word> fieldNames(1);
        fieldNames[0] = "meshDisp";

        TMToNeu.add("sourceFields", fieldNames);
        TMToNeu.add("targetFields", tgtNames);
        TMToNeu.add("fieldTypes", fieldTypes);

        NeuMappingDict.add("thermoMechanicalRegion", TMToNeu);

        dictionary deformSubDict;
        dictionary deformationDict;
        deformSubDict.add("displacementField", "disp");
        deformationDict.add("neutroRegion", deformSubDict);
        multiRegionCouplingDict.add("meshDeformation", deformationDict);
    }

    mappings.add("neutroRegion", NeuMappingDict);
}

if (solveThermoMechanics)
{
    dictionary TMMappingDict;

    if (solveThermalHydraulics and solveEnergy)
    {
        dictionary THToTM;
        List<word> tgtNames(2);
        tgtNames[0] = "TStructFromTH";
        tgtNames[1] = "TFuel";
        List<word> fieldTypes(2);
        fieldTypes[0]= "scalar";
        fieldTypes[1]= "scalar";
        List<word> fieldNames(2);
        fieldNames[0] = "T.passiveStructure";
        fieldNames[1] = linkedFuel
            ? "T.cladAvForNeutronics"
            : "T.fuelAvForNeutronics";

        THToTM.add("sourceFields", fieldNames);
        THToTM.add("targetFields", tgtNames);
        THToTM.add("fieldTypes", fieldTypes);

        TMMappingDict.add("fluidRegion", THToTM);
    }

    if (solveNeutronics)
    {
        dictionary NeuToTM;
        List<word> tgtNames(1);
        tgtNames[0] = "powerDensityNeutronics";
        List<word> fieldTypes(1);
        fieldTypes[0]= "scalar";
        List<word> fieldNames(1);
        fieldNames[0] = liquidFuel
            ? "secondaryPowerDensity"
            : "powerDensity";
        NeuToTM.add("sourceFields", fieldNames);
        NeuToTM.add("targetFields", tgtNames);
        NeuToTM.add("fieldTypes", fieldTypes);
        TMMappingDict.add("neutroRegion", NeuToTM);
    }

    mappings.add("thermoMechanicalRegion", TMMappingDict);
}
multiRegionCouplingDict.add("mappings", mappings);


// Write new dictionaries
controlDict.regIOobject::write();
multiRegionCouplingDict.regIOobject::write();
