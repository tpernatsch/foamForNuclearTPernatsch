const volVectorField& lRe(FSPair_->lRe());
const volVectorField& lDh(FSPair_->structure().lDh());

scalar maxTr2PrRatio
(
	mesh_.solutionDict().subDict("PIMPLE").lookupOrDefault<scalar>
    (
        "maxTransverseToPrincipalDragRatio", 
        100.0
    )
);

forAll(cellList_, i)
{
	label celli(cellList_[i]);
    const vector& lRei(lRe[celli]);
    const vector& lDhi(lDh[celli]);
    const scalar& halfAlphaRhoMagUi((*halfAlphaRhoMagU_)[i]);
    tensor& Kdi(Kd[celli]);
    scalar maxKd(maxTr2PrRatio*Kdi[dpa_*4]);

    Kdi[dta0_*4] = min(halfAlphaRhoMagUi*0.316*pow(lRei[dta0_],-0.25)/lDhi[dta0_], maxKd);
	Kdi[dta1_*4] = min(halfAlphaRhoMagUi*0.316*pow(lRei[dta1_],-0.25)/lDhi[dta1_], maxKd);
}
