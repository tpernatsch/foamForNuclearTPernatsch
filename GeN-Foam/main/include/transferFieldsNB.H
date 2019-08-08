

TNB.primitiveFieldRef() = thermo.T().primitiveField();
TNB.correctBoundaryConditions();

TavFuelNB.primitiveFieldRef() = TavFuel.primitiveField();
TavFuelNB.correctBoundaryConditions();

TavCladNB.primitiveFieldRef() = TavClad.primitiveField();
TavCladNB.correctBoundaryConditions();

TstructuresNB.primitiveFieldRef() = Tstructures.primitiveField();
TstructuresNB.correctBoundaryConditions();

porosityNB.primitiveFieldRef() = porousMedium.voidFraction().primitiveField();
porosityNB.correctBoundaryConditions();

rhoNB.primitiveFieldRef() = rho.primitiveField();
rhoNB.correctBoundaryConditions();

rhokNB.primitiveFieldRef() = rhok.primitiveField();
rhokNB.correctBoundaryConditions();

UNB.primitiveFieldRef() = U.primitiveField();
UNB.correctBoundaryConditions();

muNB.primitiveFieldRef() = thermo.mu()().primitiveField();
muNB.correctBoundaryConditions();

alphatNB.primitiveFieldRef() = turb.alphat()().primitiveField();
alphatNB.correctBoundaryConditions();
