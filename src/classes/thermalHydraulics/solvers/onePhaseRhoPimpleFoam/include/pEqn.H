if (!pimple_.SIMPLErho())
{
    rho_ = thermo_.rho();
}

// Thermodynamic density needs to be updated by thermo_.psi()*d(p_) after the
// pressure solution
const volScalarField psip0(thermo_.psi()*p_);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho_*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_));

if (pimple_.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho_)*fvc::flux(HbyA)
  + MRF_.zeroFilter(rhorAUf*fvc::ddtCorr(rho_, U_, phi_, rhoUf_))
);

fvc::makeRelative(phiHbyA, rho_, U_);
MRF_.makeRelative(fvc::interpolate(rho_), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_, rho_, U_, phiHbyA, rhorAUf, MRF_);

if (pimple_.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(thermo_.psi())/fvc::interpolate(rho_))*phiHbyA
    );

    phiHbyA -= fvc::interpolate(thermo_.psi()*p_)*phiHbyA/fvc::interpolate(rho_);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho_) + thermo_.psi()*correction(fvm::ddt(p_))
      + fvc::div(phiHbyA) + fvm::div(phid, p_)
     ==
        fvOptions_(thermo_.psi(), p_, rho_.name())
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p_));

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.solve(p_.select(pimple_.finalInnerIter()));

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho_) + thermo_.psi()*correction(fvm::ddt(p_))
      + fvc::div(phiHbyA)
     ==
        fvOptions_(thermo_.psi(), p_, rho_.name())
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAUf, p_));

        pEqn.solve(p_.select(pimple_.finalInnerIter()));

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqn.flux();
        }
    }
}


// Explicitly relax pressure for momentum corrector
p_.relax();

U_ = HbyA - rAU*fvc::grad(p_);
U_.correctBoundaryConditions();
fvOptions_.correct(U_);
K_ = 0.5*magSqr(U_);

if (pressureControl_.limit(p_))
{
    p_.correctBoundaryConditions();
}

thermo_.correctRho(thermo_.psi()*p_ - psip0, rhoMin_, rhoMax_) ;

#include "rhoEqn.H"
compressibleContinuityErrs();

rho_ = thermo_.rho();

// Correct rhoUf_ if the mesh_ is moving
fvc::correctRhoUf(rhoUf_, rho_, U_, phi_);

if (thermo_.dpdt())
{
    dpdt_ = fvc::ddt(p_);

    if (mesh_.moving())
    {
        dpdt_ -= fvc::div(fvc::meshPhi(rho_, U_), p_);
    }
}
