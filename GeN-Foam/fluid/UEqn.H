    // Solve the Momentum equation
    

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(rho, U)
      + (1/porosity)*fvm::div(phi, U)
	  + turb.divDevRhoReff(U)
	  - porousMedium.semiImplicitMomentumSource(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();
    UEqn.relax();

    //"reconstruct" stuff is used to fake a staggered grid (rhie chow)
    if (momentumPredictor)
    {
		if(solveFluidMechanics)
		{
			Info<< "Solve the Momentum equation " << nl << endl;
			solve
			(
				UEqn
			 	==
				porosity*fvc::reconstruct
				(
					(
					  - ghf*fvc::snGrad(rho*rhok)//p_rgh=p-rho*rhok*g*z
					  - fvc::snGrad(p_rgh)
					)*fluidMesh.magSf()
				),
				fluidMesh.solver(U.select(finalIter))

			);

			realU = U / porosity;
			realU.correctBoundaryConditions();
			K = 0.5*magSqr(realU);
		}
    }




