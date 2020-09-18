surfaceScalarField alphaf(fvc::interpolate(alpha));

//- Update fluxes
fluid_.alphaPhi() = alphaf*fluid_.phi();
fluid_.alphaRhoPhi() = 
    fvc::interpolate(rho)*fluid_.alphaPhi();

//- Check if compressible
dimensionedScalar compressibility(fvc::domainIntegrate(fluid_.thermo().psi()));
dimensionedScalar VTot("", dimVol, gSum(mesh_.V()));
bool compressible(compressibility.value()/VTot.value() >= 1e-8);

//- Construct 1/A field and interpolate to faces
volScalarField rAU
(
    "rAU",
    1.0/(UEqn.A())
);
surfaceScalarField alpharAUf
(
    fvc::interpolate(alpha*rAU)
);

//- Pressure corrector loop
while (pimple_.correct())
{
    this->correctContErr();

    p_rgh_ = p_ - rho*gh_;

    //- Correct fixed-flux BCs to be consistent with the velocity BCs
    MRF_.correctBoundaryFlux(U, fluid_.phi());

    //- Construct H/A field
    volVectorField HbyA
    (
        "HbyA",
        rAU*UEqn.H()
    );

    //- Construct flux due to explicit sources, namely buoyancy/gravity effects
    //  and the structure momentum source, if present
    //- The first term here below is the gravity term. In case one wanted to 
    //  implement a Boussinesq approximation, it should be enough to create
    //  an additional non-dimensional field (normally called rhok) and 
    //  equal to (1-Beta*DeltaT). One can then: 1) set the density as constant 
    //  in the thermoPhysicalProperties; 2) multiply the rho by rhok in the 
    //  term here below; 3) do not forget to pass to neutronics rho*rhok 
    surfaceScalarField phiSource
    (
        -alpharAUf*ghf_*fvc::snGrad(rho)*mesh_.magSf()
    );
    if (structure_.momentumSourcePtr().valid())
    {
       phiSource += alpharAUf*fvc::flux(structure_.momentumSource());
    }
    
    /* NOTE 
        
        To add a new source term S in the momentum equation, e.g.:
        
        ddt(alpha,rho,U) + ... = + S
        
        it is sufficient to add it to phiSource as:
        
        rAUf*fvc::flux(S)

        if alpha figures in S so that it can be written as:

        S = alpha*Q

        it is better to add the following to phiSource

        alpharAUf*fvc:flux(Q)

        rather than what was said before
    */

    //- Construct phase flux
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        alphaf*
        (
            fvc::flux(HbyA)
        +   phiSource
        )
    );

    //- Pressure laplacian coefficient
    surfaceScalarField rAUf
    (
        "rAUf",
        mag(alphaf*alpharAUf)
    );

    //- Update the fixedFluxPressure BCs to ensure flux consistency. Please
    //  make sure that wherever you use this BC, the U is NOT zeroGradient,
    //  as the boundary U (and thus phi) are not yet updated at this point.
    //  You get un-physical results in those circumstances. Instead, use
    //  zeroGradient for pressure were you use zeroGradient for U
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh_.boundaryFieldRef(),
        (
            phiHbyA.boundaryField()
          - (
                alphaf.boundaryField()*fluid_.phi().boundaryField()
            )
        )/(mesh_.magSf().boundaryField()*rAUf.boundaryField())
    );

    //- Construct the compressibility parts of the pressure equation
    tmp<fvScalarMatrix> pEqnComp;
    pEqnComp =
        (
            fluid_.contErr() // +dmdt in contErr
          - fvc::Sp(fvc::ddt(alpha) + fvc::div(fluid_.alphaPhi()), rho)
        )/rho
      + (alpha*fluid_.thermo().psi()/rho)*correction(fvm::ddt(p_rgh_));

    //- Cache p prior to solve for density update
    volScalarField p_rgh_0(p_rgh_);

    //- Iterate over the pressure equation to correct for non-orthogonality
    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqnIncomp
        (
            fvc::div(phiHbyA)
          - fvm::laplacian(rAUf, p_rgh_)
        );

        fvScalarMatrix pEqn(pEqnComp() + pEqnIncomp);

        pEqn.setReference
        (
            pRefCell_, 
            compressible ? getRefCellValue(p_rgh_, pRefCell_) : pRefValue_
        );

        solve
        (
            pEqn,
            mesh_.solver(p_rgh_.select(pimple_.finalInnerIter()))
        );

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqnIncomp.flux();

            fluid_.phi() = phi_/alphaf;

            fluid_.dgdt() = -pEqnComp & p_rgh_;

            p_rgh_.relax();
            U =  
                HbyA 
            +   fvc::reconstruct
                (
                    alpharAUf*pEqnIncomp.flux()/rAUf + phiSource
                );

            U.correctBoundaryConditions();
            fvOptions_.correct(U);
            fluid_.magU() = mag(U);
        }
    }
    
    //- Update and limit the static pressure
    p_ = max(p_rgh_ + rho*gh_, pMin_);

    //- Correct p_ to account for compressibility effects in closed volumes.
    //  Do this outside on the last correct() loop as otherwise the loop itself
    //  gets as stable as the Italian government
    if 
    (
        p_rgh_.needReference() //- i.e. closed volume 
    and pimple_.corrPISO() == pimple_.nCorrPISO()
    )
    {    
        //- This commented bit was originally other compressible solvers, e.g.
        //  buoyantPimpleFoam. However, it arguably does not make much sense,
        //  as the pressure was already set in setReference, so there is no
        //  need for this correction
        /*
        if (!compressible)
        {
            p_ += dimensionedScalar
            (
                "p",
                p_.dimensions(),
                pRefValue_ - getRefCellValue(p_, pRefCell_)
            );
        }*/
        if (compressible)
        {
            dimensionedScalar deltaMass
            (
                initialFluidMass_ - fvc::domainIntegrate(rho*fluid_)
            );
            
            Info<< "Compressibile relative mass defect = " 
            << (deltaMass.value()/initialFluidMass_).value() << endl;
            p_ += deltaMass/compressibility;
        }
    }

    //- Re-limit after correction
    p_ = max(p_, pMin_);

    //- Update densities from change in p_rgh
    rho += fluid_.thermo().psi()*(p_rgh_ - p_rgh_0);

    //- Correct p_rgh for consistency with p and the updated densities
    p_rgh_ = p_ - rho*gh_;
    p_rgh_.correctBoundaryConditions();
}

dpdt_ = fvc::ddt(p_);

this->correctContErr();
volScalarField contErrRel(fluid_.contErr()/rho);

Info<< "Continuity error (avg min max) ="
    << " " << contErrRel.weightedAverage(mesh_.V()).value()
    << " " << min(contErrRel).value()
    << " " << max(contErrRel).value()
    << " 1/s"  << endl;
