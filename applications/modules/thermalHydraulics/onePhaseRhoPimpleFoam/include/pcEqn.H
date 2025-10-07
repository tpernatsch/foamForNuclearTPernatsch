if (!pimple_.SIMPLErho())
{
    rho_ = thermo_.rho();
}

// Thermodynamic density needs to be updated by thermo_.psi()*d(p_) after the
// pressure solution
const volScalarField psip0(thermo_.psi()*p_);

volScalarField rAU(1.0/UEqn.A());
volScalarField rAtU(1.0/(1.0/rAU - UEqn.H1()));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U_, p_));

if (pimple_.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    (
        fvc::interpolate(rho_)*fvc::flux(HbyA)
      + MRF_.zeroFilter
        (
            fvc::interpolate(rho_*rAU)*fvc::ddtCorr(rho_, U_, phi_, rhoUf_)
        )
    )
);

fvc::makeRelative(phiHbyA, rho_, U_);
MRF_.makeRelative(fvc::interpolate(rho_), phiHbyA);

volScalarField rhorAtU("rhorAtU", rho_*rAtU);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_, rho_, U_, phiHbyA, rhorAtU, MRF_);

if (pimple_.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(thermo_.psi())/fvc::interpolate(rho_))*phiHbyA
    );

    phiHbyA +=
        fvc::interpolate(rho_*(rAtU - rAU))*fvc::snGrad(p_)*mesh_.magSf()
      - fvc::interpolate(thermo_.psi()*p_)*phiHbyA/fvc::interpolate(rho_);

    HbyA -= (rAU - rAtU)*fvc::grad(p_);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho_) + thermo_.psi()*correction(fvm::ddt(p_))
      + fvc::div(phiHbyA) + fvm::div(phid, p_)
     ==
        fvOptions_(thermo_.psi(), p_, rho_.name())
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p_));

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
    phiHbyA += fvc::interpolate(rho_*(rAtU - rAU))*fvc::snGrad(p_)*mesh_.magSf();
    HbyA -= (rAU - rAtU)*fvc::grad(p_);

    fvScalarMatrix pDDtEqn
    (
        fvc::ddt(rho_) + thermo_.psi()*correction(fvm::ddt(p_))
      + fvc::div(phiHbyA)
     ==
        fvOptions_(thermo_.psi(), p_, rho_.name())
    );

    while (pimple_.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn(pDDtEqn - fvm::laplacian(rhorAtU, p_));

        pEqn.solve(p_.select(pimple_.finalInnerIter()));

        if (pimple_.finalNonOrthogonalIter())
        {
            phi_ = phiHbyA + pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
compressibleContinuityErrs();

// Explicitly relax pressure for momentum corrector
p_.relax();

U_ = HbyA - rAtU*fvc::grad(p_);
U_.correctBoundaryConditions();
fvOptions_.correct(U_);
K_ = 0.5*magSqr(U_);

if (pressureControl_.limit(p_))
{
    p_.correctBoundaryConditions();
}

thermo_.correctRho(thermo_.psi()*p_ - psip0, rhoMin_, rhoMax_) ;
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
