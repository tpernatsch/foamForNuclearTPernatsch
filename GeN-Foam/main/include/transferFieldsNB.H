TNB.primitiveFieldRef() = thermalHydraulics.TLiquid().primitiveField();
TNB.correctBoundaryConditions();

TavFuelNB.primitiveFieldRef() = thermalHydraulics.TavFuel().primitiveField();
TavFuelNB.correctBoundaryConditions();

TavCladNB.primitiveFieldRef() = thermalHydraulics.TavClad().primitiveField();
TavCladNB.correctBoundaryConditions();

TstructuresNB.primitiveFieldRef() = thermalHydraulics.TStructure().primitiveField();
TstructuresNB.correctBoundaryConditions();

porosityNB.primitiveFieldRef() = thermalHydraulics.alphaStructure().primitiveField();
porosityNB.correctBoundaryConditions();

rhoNB.primitiveFieldRef() = thermalHydraulics.rhoLiquid()().primitiveField();
rhoNB.correctBoundaryConditions();

UNB.primitiveFieldRef() = thermalHydraulics.ULiquid().primitiveField();
UNB.correctBoundaryConditions();
