const volScalarField& TStruct_ = mesh_.lookupObjectRef<volScalarField>("T");

surfaceScalarField faceArea(fuelOrientation_ & mesh_.Sf());

fvScalarMatrix fuelDispEqn
(
    fvm::div(faceArea,fuelDisp_, "div(fuelDisp)")
    ==
    alphaFuel_*(TFuel_-TFuelRef_)
);

vector CROrientation = -fuelOrientation_;
surfaceScalarField faceAreaCR(CROrientation & mesh_.Sf());

fvScalarMatrix CRDispEqn
(
    fvm::div(faceAreaCR,CRDisp_, "div(CRDisp)")
    ==
    -alphaCR_*(TStruct_-TCRRef_)
);

Info << "Solving fuel"<< endl;
fuelDispEqn.solve();
Info << "Solving CR"<< endl;
CRDispEqn.solve();
Info << "All solved"<< endl;

volVectorField fuelDispVector_ = (fuelDisp_+CRDisp_)*fuelOrientation_;

const volVectorField& disp_ = mesh_.lookupObjectRef<volVectorField>("D");

meshDisp_ = 
    (disp_-((disp_ & fuelOrientation_) * fuelOrientation_)) 
+   fuelDispVector_ ; //- Remove any axial variation

meshDisp_.correctBoundaryConditions();