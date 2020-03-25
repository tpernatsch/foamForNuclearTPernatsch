const scalarField& alphal(liquid->primitiveField());
const scalarField& rhol(liquid->thermo().rho().primitiveField());

if (vapour != nullptr)
{
	const scalarField& alphav(vapour->primitiveField());
	const scalarField& rhov(vapour->thermo().rho().primitiveField());
	
	bafflelessAlpha.primitiveFieldRef() = alphal+alphav;

	bafflelessRhoCool.primitiveFieldRef() = 
		(alphal*rhol+alphav*rhov)/bafflelessAlpha.primitiveField();
	
	//- Should I weigh coolant temperature by mass or volume, in two-phase
	//  scenarios? I don't know, I'll use the liquid temperature for now
	bafflelessTCool.primitiveFieldRef() = 
		liquid->thermo().T().primitiveField();
}
else
{
	bafflelessAlpha.primitiveFieldRef() = alphal;

	bafflelessRhoCool.primitiveFieldRef() = rhol;

	bafflelessTCool.primitiveFieldRef() = 
		liquid->thermo().T().primitiveField();

	//- I am not adjusting this if the vapour is present, as if that is the
	//  case, I know liquidFuel is false, and U is only actually used by the 
	//  neutronics when doing liquidFuel
	bafflelessAlphat.primitiveFieldRef() = 
		liquid->turbulence().alphat()().primitiveField();
	bafflelessAlphat.correctBoundaryConditions();

	bafflelessMu.primitiveFieldRef() = 
		liquid->thermo().mu()().primitiveField();
	bafflelessMu.correctBoundaryConditions();

	bafflelessU.primitiveFieldRef() = liquid->U().primitiveField();
	bafflelessU.correctBoundaryConditions();
}

bafflelessAlpha.correctBoundaryConditions();
bafflelessRhoCool.correctBoundaryConditions();
bafflelessTCool.correctBoundaryConditions();

bafflelessTFuelAv.primitiveFieldRef() = TFuelAv->primitiveField();
bafflelessTFuelAv.correctBoundaryConditions();
bafflelessTCladAv.primitiveFieldRef() = TCladAv->primitiveField();
bafflelessTCladAv.correctBoundaryConditions();
bafflelessTStruct.primitiveFieldRef() = TStruct.primitiveField();
bafflelessTStruct.correctBoundaryConditions();
