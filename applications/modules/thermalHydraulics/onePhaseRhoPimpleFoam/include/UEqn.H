// Solve the Momentum equation

MRF_.correctBoundaryVelocity(U_);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(rho_, U_) + fvm::div(phi_, U_)
  + MRF_.DDt(rho_, U_)
  + turbulence_->divDevRhoReff(U_)
 ==
    fvOptions_(rho_, U_)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions_.constrain(UEqn);

if (pimple_.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p_));

    fvOptions_.correct(U_);
    K_ = 0.5*magSqr(U_);
}
