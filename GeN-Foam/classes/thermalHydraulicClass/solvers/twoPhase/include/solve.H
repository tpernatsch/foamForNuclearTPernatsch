//- Update continuity errors, incluisve of mass transfer
correctContErrs();

//- Correct regimes (i.e. regime marker fields)
regimeMap_->correct();

//- Correct fluid geometry, update Dhs, iA, fracs
{
    fluid1_.dispersion() *= 0.0;
    fluid2_.dispersion() *= 0.0;
    fluid1_.Dh() *= 0.0;
    fluid2_.Dh() *= 0.0;
    iA12_ *= 0.0;
    frac1_ *= 0.0;
    frac2_ *= 0.0;

    //- Set all fields according to each regime's spatial distribution
    forAllConstIter
    (
        regimeTable,
        regimeMap_->regimes(),
        regimeIter
    )
    {
        const regime& regime(regimeIter());
        if (!regime.isCurrentlyPresent()) continue;
        regime.correctFluidGeometry
        (
            iA12_,
            frac1_,
            frac2_,
            fluid1_,
            fluid2_
        );
    }
}

//- Update dimensionless numbers (Re, etc.)
FFPair_.correct();
F1SPair_.correct();
F2SPair_.correct();

//- Correct field tables for drag and heat transfer according
//  to new regime distribution in domain
{
    //- Reset all fields relevant for inter-phase coupling
    forAllIter
    (
        volTensorFieldTable,
        Kds_,
        iter
    )
    {
        volTensorField& Kd(iter());
        Kd *= 0.0;
    }
    forAllIter
    (
        volScalarFieldTable,
        htcs_,
        iter
    )
    {
        volScalarField& htc(iter());
        htc *= 0.0;
    }

    //- Set all fields according to each regime's spatial distribution
    forAllConstIter
    (
        regimeTable,
        regimeMap_->regimes(),
        regimeIter
    )
    {
        const regime& regime(regimeIter());
        if (!regime.isCurrentlyPresent()) continue;
        regime.correctDragTable(Kds_);
        regime.correctHeatTransferTable(htcs_);
        
    }
}

if (solveFluidDynamics_)
{
    //- Solve continuity equations to compute new alphas
    #include "alphaEqns.H"

    //- Construct momentum matrices
    #include "UEqns.H"

    //- Solve pressure equation and reconstruct velocities
    if (faceMomentum_)
    {
        #include "pEqnf.H"
    }
    else
    {
        #include "pEqn.H"
    }
}
if (solveEnergy_)
{
    //- Solve fluid energy equations and update structure energy sources
    #include "EEqns.H"
}
