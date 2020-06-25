const volVectorField& lRe(FSPair_->lRe());
const volVectorField& lDh(FSPair_->structure().lDh());

forAll(cellList_, i)
{
	label celli(cellList_[i]);
    const vector& lRei(lRe[celli]);
    const vector& lDhi(lDh[celli]);
    const scalar& halfAlphaRhoMagUi((*halfAlphaRhoMagU_)[i]);
    tensor& Kdi(Kd[celli]);

    Kdi[dta0_*4] = halfAlphaRhoMagUi*0.316*pow(lRei[dta0_],-0.25)/lDhi[dta0_];
	Kdi[dta1_*4] = halfAlphaRhoMagUi*0.316*pow(lRei[dta1_],-0.25)/lDhi[dta1_];
}
