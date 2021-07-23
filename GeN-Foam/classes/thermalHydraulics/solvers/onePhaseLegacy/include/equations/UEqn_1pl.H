//- Field refs
const volScalarField& alpha(fluid_);
const volScalarField& psi(fluid_.thermo().psi());

volVectorField& UDarcy(UDarcy_());
volVectorField& U(fluid_.U());
surfaceScalarField& phiDarcy(phiDarcy_());

fvVectorMatrix UEqn(UDarcy, dimMass*dimLength/dimTime/dimTime);

//- Drag between fluid and structure rotated in global reference frame
const volTensorField& Kds(FSPair_.Kd());

//- Resetting rho to the reference value in case of incompressible 
//- Boussinesq flow. It is important, because otherwise the 
//- turbulence module would use a temperature corrected density, 
//- which we do not want.
if(incompressibleTreatment_)
{
    fluid_.thermo().rho() = fixedRho_;
    fluid_.thermo().rho().correctBoundaryConditions();
}
fluid_.turbulence().correct();

volScalarField nuEff(fluid_.turbulence().nuEff()());

UEqn = 
(
        fvm::ddt(fixedRho_, UDarcy)
    +   (1/alpha)*fvm::div(phiDarcy, UDarcy, "div(phiDarcy,UDarcy)")
    -   (1/alpha)*fvm::SuSp(fluid_.contErr(), UDarcy)
    //- Diffusion term (i.e. fluid_.turbulence().divDevRhoReff(U))
    //  but re-expressed so to: 1) either allow to use a tensor
    //  diffusivity (i.e. nuEff*structure_.tortuosity()) if the 
    //  structure is tortuous; 2) if not tortuous, still allow me to
    //  specify the div interpolation scheme (i.e. "div(RhoPhiNu,UDarcy)")
    +   (
            (structure_.tortuosityPtr().valid()) ? 
            (
            -   fvm::laplacian
                (
                    fixedRho_*nuEff*structure_.tortuosity(), 
                    UDarcy
                )
            -   fvc::div
                (
                    fixedRho_*nuEff*structure_.tortuosity()
                &   dev2(T(fvc::grad(UDarcy))), 
                    "div(RhoPhiNu,UDarcy)"
                )
            ) :
            (
            -   fvm::laplacian(fixedRho_*nuEff, UDarcy)
            -   fvc::div
                (
                    fixedRho_*nuEff*dev2(T(fvc::grad(UDarcy))), 
                    "div(RhoPhiNu,UDarcy)"
                )
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


