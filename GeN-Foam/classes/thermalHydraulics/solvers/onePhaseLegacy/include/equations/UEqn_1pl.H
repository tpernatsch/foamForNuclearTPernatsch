//- Field refs
const volScalarField& alpha(fluid_);
const volScalarField& psi(fluid_.thermo().psi());

volVectorField& UDarcy(UDarcy_());
volVectorField& U(fluid_.U());
surfaceScalarField& phiDarcy(phiDarcy_());

fvVectorMatrix UEqn(UDarcy, dimMass*dimLength/dimTime/dimTime);

//- Drag between fluid and structure rotated in global reference frame
volTensorField Kds
(
    IOobject
    (
        "Kds",
        mesh_.time().timeName(),
        mesh_
    ),
    (
        structure_.Rl2g() 
    &   (*Kds_["FSPair"])
    &   structure_.Rg2l()
    ),
    zeroGradientFvPatchScalarField::typeName
);
Kds = Kds/alpha;
Kds.correctBoundaryConditions();

/*
volScalarField magKds(mag(Kds)/sqrt(3.0));
Info<< "Kds (avg min max) = "
        << magKds.weightedAverage(mesh_.V()).value()
        << " " << min(magKds).value()
        << " " << max(magKds).value()
        << " kg/m3/s" << endl;
*/

fluid_.turbulence().correct();

volTensorField nuEff
(
    structure_.tortuosity()
   *fluid_.turbulence().nuEff()()
);

UEqn = 
(
        fvm::ddt(fixedRho_, UDarcy)
    +   (1/alpha)*fvm::div(phiDarcy, UDarcy, "div(phiDarcy,UDarcy)")
        //- The following is just a re-arrangement of div(nu*grad(U))
    +   (
            -   fvm::laplacian(fixedRho_*nuEff, UDarcy)
            -   fvc::div
                (
                    fixedRho_*nuEff&dev2(T(fvc::grad(UDarcy))), 
                    "div(RhoPhiNu,UDarcy)"
                )
        )
        //- Implicit diagonal contribution of fluid-structure drag, 
        //  explicit off-diagonal contribution
    +   fvm::Sp((1.0/3.0)*tr(Kds), UDarcy) + (dev(Kds) & UDarcy)
    ==
    (
        structure_.momentumSourcePtr().valid() ? 
        (structure_.momentumSource() + fvOptions_(fixedRho_, UDarcy)) :
        fvOptions_(fixedRho_, UDarcy)
    )
);

UEqn.relax();
fvOptions_.constrain(UEqn);

solve
(
    UEqn
  ==
    alpha*fvc::reconstruct
    (
        (
          - ghf_*fvc::snGrad(fixedRho_*rhok_)
          - fvc::snGrad(p_rgh_)
        )*mesh_.magSf()
    ),
    mesh_.solver(UDarcy.select(pimple_.finalInnerIter()))
);

U = UDarcy / alpha;
U.correctBoundaryConditions();

fvOptions_.correct(UDarcy);


