

Info<< "Reading thermophysical properties \n" << endl;


autoPtr<rhoThermo> pThermo
(
    rhoThermo::New(fluidMesh)
);
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

Info<< "Done \n" << endl;


Info<< "Reading fluid fields \n" << endl;

uniformDimensionedScalarField beta
(
	IOobject
	(
		"beta",
		runTime.constant(),
		fluidMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

uniformDimensionedScalarField Tref
(
	IOobject
	(
		"Tref",
		runTime.constant(),
		fluidMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

uniformDimensionedScalarField rhoRef
(
	IOobject
	(
		"rhoRef",
		runTime.constant(),
		fluidMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

uniformDimensionedVectorField g
(
	IOobject
	(
		"g",
		runTime.constant(),
		fluidMesh,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

volScalarField rho
(
	IOobject
	(
		"rho",
		runTime.timeName(),
		fluidMesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	thermo.rho()
);

const volScalarField& psi = thermo.psi();


volScalarField rhok
(
	IOobject
	(
		"rhok",
		runTime.timeName(),
		fluidMesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	fluidMesh
);

volVectorField U
(
	IOobject
	(
		"U",
		runTime.timeName(),
		fluidMesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	fluidMesh
);

volVectorField realU
(
	IOobject
	(
		"realU",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	U
);

surfaceScalarField phi
(
	IOobject
	(
		"phi",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	linearInterpolate(rho*U) & fluidMesh.Sf()
);

volScalarField p_rgh
(
	IOobject
	(
		"p_rgh",
		runTime.timeName(),
		fluidMesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	fluidMesh
);

volScalarField gh("gh",	g & fluidMesh.C());

p_rgh = thermo.p() - rho*gh; // Force p_rgh to be consistent with p

surfaceScalarField ghf("ghf", g & fluidMesh.Cf());

volScalarField K("K", 0.5*magSqr(realU));

volScalarField dpdt
(
	IOobject
	(
		"dpdt",
		runTime.timeName(),
		fluidMesh
	),
	fluidMesh,
	dimensionedScalar
	(
		"dpdt",
		thermo.p().dimensions()/dimTime,
		0
	)	
);

volScalarField Re
(
	IOobject
	(
		"Re",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        fluidMesh,
	dimensionedScalar("", dimless, 0.0),
	zeroGradientFvPatchScalarField::typeName
);

volScalarField volFuelPower
(
	IOobject
	(
		"volFuelPower",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        fluidMesh,
	dimensionedScalar("", dimensionSet(1,-1,-3,0,0,0,0), 0.0),
	zeroGradientFvPatchScalarField::typeName
);

volScalarField TavFuel
(
	IOobject
	(
		"TavFuel",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        thermo.T()
);

volScalarField TavClad
(
	IOobject
	(
		"TavClad",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        thermo.T()
);

volScalarField TinnerFuel
(
	IOobject
	(
		"TinnerFuel",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        thermo.T()
);

volScalarField TouterFuel
(
	IOobject
	(
		"TouterFuel",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        thermo.T()
);

volScalarField TinnerClad
(
	IOobject
	(
		"TinnerClad",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        thermo.T()
);

volScalarField TouterClad
(
	IOobject
	(
		"TouterClad",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        thermo.T()
);

volScalarField Tstructures
(
	IOobject
	(
		"Tstructures",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
        thermo.T()
);

Info<< "Done \n" << endl;

Info<< "Reading main fluid classes \n" << endl;

autoPtr<porousMedium> pPorousMedium
(
	porousMedium::New
	(
		rho,
		realU,
		thermo
	)
);
porousMedium& porousMedium = pPorousMedium();

autoPtr<subscaleFuel> pSubscaleFuel
(
	subscaleFuel::New
	(
		fluidMesh
	)
);
subscaleFuel& subscaleFuel = pSubscaleFuel();


autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
compressible::turbulenceModel& turb = turbulence();

//displacement field for the mesh (should come from mechanical solver)
volVectorField DFluid
(
	IOobject
	(
		"DFluid",
		runTime.timeName(),
		fluidMesh,
		IOobject::READ_IF_PRESENT,
		IOobject::AUTO_WRITE
	),
	fluidMesh,
	dimensionedVector("d_zero", dimensionSet(0,1,0,0,0,0,0), vector(0,0,0)),
	zeroGradientFvPatchScalarField::typeName
);

Info<< "Done \n" << endl;

