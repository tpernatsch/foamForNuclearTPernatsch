
//- Pressure corrector loop
while (pimple_.correct())
{
    //- Check if compressible
    dimensionedScalar compressibility(fvc::domainIntegrate(fluid_.thermo().psi()));
    dimensionedScalar VTot("", dimVol, gSum(mesh_.V()));
    bool compressible(compressibility.value()/VTot.value() >= 1e-8);
    
    //- Construct 1/A field and interpolate to faces
    volScalarField rAU("rAUDarcy", 1.0/UEqn.A());
    surfaceScalarField alphaRhorAUf(fvc::interpolate(alpha*rho*rAU));
    
    volVectorField HbyA("HbyA", UDarcy);
    HbyA = rAU*(UEqn.H());
    
    surfaceScalarField phig(-alphaRhorAUf*ghf_*fvc::snGrad(rho*rhok)*mesh_.magSf());
    
    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        (
            fvc::flux(rho*HbyA)
            + alphaRhorAUf*fvc::ddtCorr(rho, UDarcy, phiDarcy)
        )
        + phig
    );

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, UDarcy, phiHbyA, alphaRhorAUf);

    // fvScalarMatrix p_rghDDtEqn
    // (
    //     // (compr)*(fvc::ddt(porosity,rho) + psi*correction(fvm::ddt(porosity,p_rgh)))
    //     + fvc::div(phiHbyA)
    // );

    //- Iterate over the pressure equation to correct for non-orthogonality
    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            /*p_rghDDtEqn*/
            - fvm::laplacian(alphaRhorAUf, p_rgh_)
            + fvc::div(phiHbyA)
        );

        if(!compressible)
        {
            pEqn.setReference
            (
                pRefCell_, 
                compressible ? getRefCellValue(p_rgh_, pRefCell_) : pRefValue_
            );
        }
        
        solve
        (
            pEqn,
            mesh_.solver(p_rgh_.select(pimple_.finalInnerIter()))
        );

        if (pimple_.finalNonOrthogonalIter())
        {
            phiDarcy = phiHbyA + pEqn.flux();

            UDarcy = HbyA + rAU*alpha*fvc::reconstruct((phig + pEqn.flux())/alphaRhorAUf);

            UDarcy.correctBoundaryConditions();
            
            U = UDarcy/alpha;

            U.correctBoundaryConditions();

            fvOptions_.correct(UDarcy);

            fluid_.magU() = mag(U);
        }
    }
    
    //- Update and limit the static pressure
    p_ = max(p_rgh_ + rho*rhok*gh_, pMin_);

    // if 
    // (
    //     p_rgh_.needReference() //- i.e. closed volume 
    // and pimple_.corrPISO() == pimple_.nCorrPISO()
    // )
    // {    
    //     if (compressible)
    //     {
    //         dimensionedScalar deltaMass
    //         (
    //             initialFluidMass_ - fvc::domainIntegrate(rho*fluid_)
    //         );
    //         
    //         Info<< "Compressibile relative mass defect = " 
    //         << (deltaMass.value()/initialFluidMass_).value() << endl;
    //         p_ += deltaMass/compressibility;
    //     }
    // }
// 
    // //- Re-limit after correction
    // p_ = max(p_, pMin_);
// 
    // //- Update densities from change in p_rgh
    // rho += fluid_.thermo().psi()*(p_rgh_ - p_rgh_0);
// 
    // //- Correct p_rgh for consistency with p and the updated densities
    // p_rgh_ = p_ - rho*gh_;
    // p_rgh_.correctBoundaryConditions();
}

dpdt_ = fvc::ddt(p_);

this->correctContErr();

scalar magContErr = runTime_.deltaTValue()*
        mag(fluid_.contErr())().weightedAverage(mesh_.V()).value();

scalar contErr = runTime_.deltaTValue()*
        fluid_.contErr().weightedAverage(mesh_.V()).value();

cumulativeContErr_ += magContErr;

Info<< "Continuity error ="
    << " Int. Mag.: " << magContErr
    << " Int.: " << contErr
    << " Cumulative: " << cumulativeContErr_
    << endl;
