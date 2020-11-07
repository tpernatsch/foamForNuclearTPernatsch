//- The only thing that should be defined by the class that includes this is
//  the scalar drag factor fd(Re)

//- Alpha is the void fraction
volScalarField alpha(1.0-FSPair_->structure());
const volScalarField& magU(FSPair_->fluidRef().magU());
const volScalarField& rho(FSPair_->fluidRef().thermo().rho());

//- Minimum allowable velocity in the calculation Kd to avoid 0 Kd
scalar minMagU
(
    mesh_.solutionDict().subDict("PIMPLE").lookupOrDefault<scalar>
    (
        "minMagU", 
        0.0
    )
);

if (isotropic_) //- Use scalar values for Re, Dh, compute a single Kd and
                //  apply it to all directions
{
    const volScalarField& Re(FSPair_->Re());
    const volScalarField& Dh(FSPair_->structure().Dh());
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        const scalar& Rei(Re[celli]);
        const scalar& Dhi(Dh[celli]);

        tensor& Kdi(Kd[celli]);
        
        //- Yes, i, not celli, as alphaRhoMagU_ is not a volScalarField
        scalar& halfAlphaRhoMagUi((*halfAlphaRhoMagU_)[i]);
        halfAlphaRhoMagUi = 
            (0.5*alpha[celli]*rho[celli]*max(magU[celli], minMagU));
        
        //- Same fd, same value for all components
        scalar KdValue = (halfAlphaRhoMagUi/Dhi)*fd(Rei);
        
        Kdi[dpa_*4] = KdValue;
        Kdi[dta0_*4] = KdValue;
        Kdi[dta1_*4] = KdValue;
    }
}
else if (transverseDragModel_ == transverseDragModel::same)
                //- Use vector values for Re, Dh for each direction and
                //  compute a scalar Kd value for each component, each
                //  obtained via the fd correction
{
    const volVectorField& lRe(FSPair_->lRe());
    const volVectorField& lDh(FSPair_->structure().lDh());
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        const vector& lRei(lRe[celli]);
        const vector& lDhi(lDh[celli]);
        tensor& Kdi(Kd[celli]);
        
        //- Yes, i, not celli, as alphaRhoMagU_ is not a volScalarField
        scalar& halfAlphaRhoMagUi((*halfAlphaRhoMagU_)[i]);
        halfAlphaRhoMagUi = 
            (0.5*alpha[celli]*rho[celli]*max(magU[celli], minMagU));
        
        //- Rehme correlation for all direction with directional Re, Dh
        Kdi[dpa_*4] = 
            (halfAlphaRhoMagUi/lDhi[dpa_])*fd(lRei[dpa_]);
        Kdi[dta0_*4] = 
            (halfAlphaRhoMagUi/lDhi[dta0_])*fd(lRei[dta0_]);
        Kdi[dta1_*4] = 
            (halfAlphaRhoMagUi/lDhi[dta1_])*fd(lRei[dta1_]);
    }
}
else            //- Use vector values for Re, Dh for each direction and
                //  compute a scalar Kd value for each component, the
                //  component in the principal axis is obtained via the
                //  fd correlation and the other transverse components
                //  are computed via the selected transverseDragModel
{
    const volVectorField& lRe(FSPair_->lRe());
    const volVectorField& lDh(FSPair_->structure().lDh());
    forAll(cellList_, i)
    {
        label celli(cellList_[i]);
        const vector& lRei(lRe[celli]);
        const vector& lDhi(lDh[celli]);
        tensor& Kdi(Kd[celli]);
        
        //- Yes, i, not celli, as alphaRhoMagU_ is not a volScalarField
        scalar& halfAlphaRhoMagUi((*halfAlphaRhoMagU_)[i]);
        halfAlphaRhoMagUi = 
            (0.5*alpha[celli]*rho[celli]*max(magU[celli], minMagU));
        
        //- Rehme correlation in the principal direction
        Kdi[dpa_*4] = 
            (halfAlphaRhoMagUi/lDhi[dpa_])*fd(lRei[dpa_]);
    }
    correctTransverseKd(Kd);
}
