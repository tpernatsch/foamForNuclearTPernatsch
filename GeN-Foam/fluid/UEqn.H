// Solve the Momentum equation


tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho, U)
  + (1/porosity)*fvm::div(phi, U)
  + turb.divDevRhoReff(U)
  - porousMedium.semiImplicitMomentumSource(U)
 ==
    fvOptions(rho, U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

//"reconstruct" used to fake a staggered grid (rhie chow)
if (momentumPredictor)
{
    if(solveFluidMechanics)
    {
        Info<< "Solve the Momentum equation " << nl << endl;
        
        solvPerfU =
            solve
            (
                UEqn
              ==
                porosity*fvc::reconstruct
                (
                    (
                      - ghf*fvc::snGrad(rho*rhok)
                      - fvc::snGrad(p_rgh)
                    )*fluidMesh.magSf()
                ),
                fluidMesh.solver(U.select(finalIter))
            );

        realU = U / porosity;
        realU.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(realU);
    }
}




