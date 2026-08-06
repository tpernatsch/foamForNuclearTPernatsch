MRF_.correctBoundaryVelocity(U_);

fvVectorMatrix UEqn
(
    fvm::ddt(U_) + fvm::div(phi_, U_)
  + MRF_.DDt(U_)
  + turbulence_->divDevReff(U_)
 ==
    FP_.SU(U_, true)
  + fvOptions_(U_)
);

UEqn.relax();

fvOptions_.constrain(UEqn);

if (pimple_.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
              - ghf_*fvc::snGrad(rhok_)
              - fvc::snGrad(p_rgh_)
            )*mesh_.magSf()
        )
    );

    fvOptions_.correct(U_);
}
