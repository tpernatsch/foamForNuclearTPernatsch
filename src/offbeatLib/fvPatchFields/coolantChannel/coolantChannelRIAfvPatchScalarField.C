/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "coolantChannelRIAfvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IF97.H"
#include "globalOptions.H"
#include "wedgePolyPatch.H"
#include "mathematicalConstants.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * //

scalarField coolantChannelRIAfvPatchScalarField::calculateFaceHeights() const
{
    const fvPatch& p = this->patch();
    const pointField& localPoints = p.patch().localPoints();

    const vector pinDir(
        p.patch().boundaryMesh().mesh().foundObject<globalOptions>("globalOptions")?
        p.patch().boundaryMesh().mesh().lookupObject<globalOptions>("globalOptions").pinDirection() :
        vector(0, 0, 1) // Default pin direction if globalOptions or pinDirection is not found
    );

    scalarField faceHeights(p.size(), 0.0);

    forAll(p.patch().localFaces(), faceI)
    {
        const labelList& pointLabels = p.patch().localFaces()[faceI];

        scalar zMin = GREAT;
        scalar zMax = -GREAT;

        forAll(pointLabels, pointI)
        {
            const scalar z = localPoints[pointLabels[pointI]] & pinDir;
            zMin = min(zMin, z);
            zMax = max(zMax, z);
        }

        faceHeights[faceI] = zMax - zMin;
    }

    return faceHeights;
}


void coolantChannelRIAfvPatchScalarField::validateGeometry() const
{
    const polyBoundaryMesh& patches =
        patch().boundaryMesh().mesh().boundaryMesh();

    bool hasWedge = false;
    forAll(patches, patchI)
    {
        if (patches[patchI].size() && isType<wedgePolyPatch>(patches[patchI]))
        {
            hasWedge = true;
            break;
        }
    }

    if (!hasWedge)
    {
        FatalErrorInFunction
            << "The coolantChannelRIA BC requires wedge patches so that the "
            << "channel surface area can be interpreted consistently. "
            << "No wedge patch was found."
            << abort(FatalError);
    }

    if (coolant_.A() <= SMALL)
    {
        FatalErrorInFunction
            << "flowArea must be > 0. Current value = " << coolant_.A()
            << abort(FatalError);
    }

    if (coolant_.DH() <= SMALL)
    {
        FatalErrorInFunction
            << "hydraulicDiameter must be > 0. Current value = "
            << coolant_.DH()
            << abort(FatalError);
    }
}


void coolantChannelRIAfvPatchScalarField::updateFaceHeightsIfNeeded()
{
    if
    (
        this->db().time().value() > timeCheck_
     && db().foundObject<fvMesh>("referenceMesh")
    )
    {
        faceHeights_ = calculateFaceHeights();
    }
}


scalarField coolantChannelRIAfvPatchScalarField::computeSolidSideConductance()
const
{
    const fvPatch& p = this->patch();
    const scalarField& deltaCoeffs = p.deltaCoeffs();

    scalarField alphaP(this->size(), 0.0);

    if (this->patch().boundaryMesh().mesh().foundObject<volTensorField>(kappaName_))
    {
        const fvPatchField<tensor>& k =
            p.lookupPatchField<volTensorField, tensor>(kappaName_);

        const vectorField n = p.Sf()/p.magSf();
        alphaP = ((n & k) & n) * deltaCoeffs;
    }
    else
    {
        const fvPatchField<scalar>& k =
            p.lookupPatchField<volScalarField, scalar>(kappaName_);

        alphaP = k * deltaCoeffs;
    }

    return alphaP;
}


void coolantChannelRIAfvPatchScalarField::computeHeatFlux
(
    const scalarField& alphaP
)
{
    heatFlux_ = (patchInternalField() - *this) * alphaP;
}


void coolantChannelRIAfvPatchScalarField::advanceTimeStep()
{
    if (this->db().time().value() <= timeCheck_)
    {
        return;
    }

    vaporizedThicknessOld_ = vaporizedThickness_;
    conductionTimeOld_ = conductionTime_;
    filmBoilingTimeOld_ = filmBoilingTime_;
    timeCheck_ = this->db().time().value();
}


scalar coolantChannelRIAfvPatchScalarField::computeVaporizedThickness
(
    const label faceI,
    const scalar heatFlux,
    const scalar hSP,
    const scalar Tw,
    const scalar T0,
    const scalar p,
    const scalar deltaT
) const
{
    const scalar Tfilm = 0.5*(Tw + T0);
    const scalar netFlux = heatFlux - hSP*(Tw - T0);
    const scalar Hlv = IF97::latentHeatVaporization_p(p);
    const scalar rho = IF97::rhomass_Tp(Tfilm, p);

    return max
    (
        0.0,
        vaporizedThicknessOld_[faceI]
      + (netFlux / max(Hlv, SMALL) / max(rho, SMALL)) * deltaT
    );
}


scalar coolantChannelRIAfvPatchScalarField::computeTransitionDecay
(
    const scalar qCHF,
    const scalar qFilmLeid,
    const scalar qFilmDeriv,
    const scalar lDT_CHF
) const
{
    // Derive Bessiron decay constant k from the minimum condition
    // dΦ/dTw = 0 at TLeidenfrost, where Φ = qFloor + (qCHF-qFloor)*exp(-k*x²)
    // and x = Tw - TCHF.  Setting dΦ/dTw = 0 at x = lDT_CHF gives:
    //   qFilmDeriv = 2*k*lDT_CHF*(qCHF - qFloor)*exp(-k*lDT_CHF²)
    // Substituting u = k*lDT_CHF² and C = qFilmDeriv*lDT_CHF/(2*(qCHF-qFloor)):
    //   g(u) = C*(exp(u) - 1) - u = 0
    const scalar dq = qCHF - qFilmLeid;
    if (dq < SMALL || qFilmDeriv < SMALL) return SMALL;

    const scalar C = qFilmDeriv*lDT_CHF / (2.0*dq);
    if (C >= 1.0) return SMALL;

    // Newton-Raphson on g(u) = C*(exp(u)-1) - u = 0
    // Initial guess from large-u asymptotic: u ~ log(u/C) ~ 2*log(1/C)
    scalar u = max(2.0*log(1.0/max(C, SMALL)), 1.0);
    for (label iter = 0; iter < 50; ++iter)
    {
        const scalar eu = exp(min(u, 700.0));
        const scalar g  = C*(eu - 1.0) - u;
        const scalar gp = C*eu - 1.0;
        if (mag(gp) < SMALL) break;
        const scalar du = -g/gp;
        u += du;
        u  = max(u, SMALL);
        if (mag(du) < 1e-12*u) break;
    }

    return max(u/sqr(lDT_CHF), SMALL);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

coolantChannelRIAfvPatchScalarField::coolantChannelRIAfvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    kappaName_("k"),
    correlations_(),
    singlePhase_("ChurchillChu"),
    filmBoiling_("Sakurai"),
    coolant_(p.size()),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    relax_(1.0),
    timeCheck_(0.0),
    heatFlux_(p.size(), 0.0),
    faceHeights_(p.size(), 0.0),
    filmBoilingRIACorrection_(5.0),
    filmBoilingRIACorrectionDuration_(15.0),
    leidenfrostdTBessiron_(450.0),
    onbdTBessiron_(20.0),
    thicknessBessironRIA_(3e-5),
    resetDeltaTBessiron_(1.0),
    pulseStartTime_(0.0),
    vaporizedThickness_(p.size(), 0.0),
    vaporizedThicknessOld_(p.size(), 0.0),
    qCHF_(p.size(), 0.0),
    qLeidenfrost_(p.size(), 0.0),
    kBessiron_(p.size(), 0.0),
    conductionTime_(p.size(), 0.0),
    conductionTimeOld_(p.size(), 0.0),
    filmBoilingTime_(p.size(), 0.0),
    filmBoilingTimeOld_(p.size(), 0.0)
{}


coolantChannelRIAfvPatchScalarField::coolantChannelRIAfvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict,
    const bool valueRequired
)
:
    fvPatchField<scalar>(p, iF, dict, valueRequired),
    kappaName_(dict.lookupOrDefault<word>("kappa", "k")),
    correlations_(dict.subDict("correlations")),
    singlePhase_
    (
        correlations_.lookupOrDefault<word>("singlePhase", "ChurchillChu")
    ),
    filmBoiling_
    (
        correlations_.lookupOrDefault<word>("filmBoiling", "Sakurai")
    ),
    coolant_(dict, p.size(), false, false, false),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    relax_(dict.lookupOrDefault<scalar>("relax", 1)),
    timeCheck_(0.0),
    heatFlux_(p.size(), 0.0),
    faceHeights_(p.size(), 0.0),
    filmBoilingRIACorrection_
    (
        dict.lookupOrDefault<scalar>("filmBoilingRIACorrection", 5.0)
    ),
    filmBoilingRIACorrectionDuration_
    (
        dict.lookupOrDefault<scalar>("filmBoilingRIACorrectionDuration", 15.0)
    ),
    leidenfrostdTBessiron_
    (
        dict.lookupOrDefault<scalar>("leidenfrostdTBessiron", 450.0)
    ),
    onbdTBessiron_
    (
        dict.lookupOrDefault<scalar>("onbdTBessiron", 20.0)
    ),
    thicknessBessironRIA_
    (
        dict.lookupOrDefault<scalar>("thicknessBessironRIA", 3e-5)
    ),
    resetDeltaTBessiron_
    (
        dict.lookupOrDefault<scalar>("resetDeltaTBessiron", 1.0)
    ),
    pulseStartTime_
    (
        dict.lookupOrDefault<scalar>("pulseStartTime", 0.0)
    ),
    vaporizedThickness_
    (
        dict.found("vaporizedThickness")
      ? scalarField("vaporizedThickness", dict, p.size())
      : scalarField(p.size(), 0.0)
    ),
    vaporizedThicknessOld_(p.size(), 0.0),
    qCHF_
    (
        dict.found("qCHF")
      ? scalarField("qCHF", dict, p.size())
      : scalarField(p.size(), 0.0)
    ),
    qLeidenfrost_
    (
        dict.found("qLeidenfrost")
      ? scalarField("qLeidenfrost", dict, p.size())
      : scalarField(p.size(), 0.0)
    ),
    kBessiron_
    (
        dict.found("kBessiron")
      ? scalarField("kBessiron", dict, p.size())
      : scalarField(p.size(), 0.0)
    ),
    conductionTime_
    (
        dict.found("conductionTime")
      ? scalarField("conductionTime", dict, p.size())
      : scalarField(p.size(), 0.0)
    ),
    conductionTimeOld_(p.size(), 0.0),
    filmBoilingTime_
    (
        dict.found("filmBoilingTime")
      ? scalarField("filmBoilingTime", dict, p.size())
      : scalarField(p.size(), 0.0)
    ),
    filmBoilingTimeOld_(p.size(), 0.0)
{
    // if (coolant_.enthalpyModel())
    // {
    //     FatalErrorInFunction
    //         << "coolantChannelRIA does not support enthalpyModel. "
    //         << "Provide an axialProfileDict instead."
    //         << abort(FatalError);
    // }

    const scalar userTime =
        patch().boundaryMesh().mesh().time().timeToUserTime
        (
            this->db().time().value()
        );

    coolant_.updateBoundaryInputs(userTime);
    faceHeights_ = calculateFaceHeights();
    validateGeometry();

    vaporizedThicknessOld_ = vaporizedThickness_;
    conductionTimeOld_ = conductionTime_;
    filmBoilingTimeOld_ = filmBoilingTime_;
}


coolantChannelRIAfvPatchScalarField::coolantChannelRIAfvPatchScalarField
(
    const coolantChannelRIAfvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<scalar>(ptf, p, iF, mapper),
    kappaName_(ptf.kappaName_),
    correlations_(ptf.correlations_),
    singlePhase_(ptf.singlePhase_),
    filmBoiling_(ptf.filmBoiling_),
    coolant_(ptf.coolant_),
    h_(ptf.h_),
    alpha_(ptf.alpha_),
    relax_(ptf.relax_),
    timeCheck_(ptf.timeCheck_),
    heatFlux_(ptf.heatFlux_),
    faceHeights_(ptf.faceHeights_),
    filmBoilingRIACorrection_(ptf.filmBoilingRIACorrection_),
    filmBoilingRIACorrectionDuration_(ptf.filmBoilingRIACorrectionDuration_),
    leidenfrostdTBessiron_(ptf.leidenfrostdTBessiron_),
    onbdTBessiron_(ptf.onbdTBessiron_),
    thicknessBessironRIA_(ptf.thicknessBessironRIA_),
    resetDeltaTBessiron_(ptf.resetDeltaTBessiron_),
    pulseStartTime_(ptf.pulseStartTime_),
    vaporizedThickness_(ptf.vaporizedThickness_),
    vaporizedThicknessOld_(ptf.vaporizedThicknessOld_),
    qCHF_(ptf.qCHF_),
    qLeidenfrost_(ptf.qLeidenfrost_),
    kBessiron_(ptf.kBessiron_),
    conductionTime_(ptf.conductionTime_),
    conductionTimeOld_(ptf.conductionTimeOld_),
    filmBoilingTime_(ptf.filmBoilingTime_),
    filmBoilingTimeOld_(ptf.filmBoilingTimeOld_)
{}


coolantChannelRIAfvPatchScalarField::coolantChannelRIAfvPatchScalarField
(
    const coolantChannelRIAfvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(ptf, iF),
    kappaName_(ptf.kappaName_),
    correlations_(ptf.correlations_),
    singlePhase_(ptf.singlePhase_),
    filmBoiling_(ptf.filmBoiling_),
    coolant_(ptf.coolant_),
    h_(ptf.h_),
    alpha_(ptf.alpha_),
    relax_(ptf.relax_),
    timeCheck_(ptf.timeCheck_),
    heatFlux_(ptf.heatFlux_),
    faceHeights_(ptf.faceHeights_),
    filmBoilingRIACorrection_(ptf.filmBoilingRIACorrection_),
    filmBoilingRIACorrectionDuration_(ptf.filmBoilingRIACorrectionDuration_),
    leidenfrostdTBessiron_(ptf.leidenfrostdTBessiron_),
    onbdTBessiron_(ptf.onbdTBessiron_),
    thicknessBessironRIA_(ptf.thicknessBessironRIA_),
    resetDeltaTBessiron_(ptf.resetDeltaTBessiron_),
    pulseStartTime_(ptf.pulseStartTime_),
    vaporizedThickness_(ptf.vaporizedThickness_),
    vaporizedThicknessOld_(ptf.vaporizedThicknessOld_),
    qCHF_(ptf.qCHF_),
    qLeidenfrost_(ptf.qLeidenfrost_),
    kBessiron_(ptf.kBessiron_),
    conductionTime_(ptf.conductionTime_),
    conductionTimeOld_(ptf.conductionTimeOld_),
    filmBoilingTime_(ptf.filmBoilingTime_),
    filmBoilingTimeOld_(ptf.filmBoilingTimeOld_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coolantChannelRIAfvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    
    updateFaceHeightsIfNeeded();

    const scalar userTime =
        patch().boundaryMesh().mesh().time().timeToUserTime
        (
            this->db().time().value()
        );

    coolant_.updateBoundaryInputs(userTime);

    const fvPatch& p = this->patch();
    const vector pinDir(
        p.patch().boundaryMesh().mesh().foundObject<globalOptions>("globalOptions")?
        p.patch().boundaryMesh().mesh().lookupObject<globalOptions>("globalOptions").pinDirection() :
        vector(0, 0, 1) // Default pin direction if globalOptions or pinDirection is not found
    );

    const scalarField zetas = patch().Cf() & pinDir;
    const scalarField radii = mag(patch().Cf() - zetas*pinDir);

    const scalarField alphaP = computeSolidSideConductance();
    computeHeatFlux(alphaP);

    const scalar deltaT =
        patch().boundaryMesh().mesh().time().deltaTValue();

    coolant_.updateBulkState
    (
        userTime,
        zetas,
        radii,
        faceHeights_,
        heatFlux_,
        deltaT
    );

    // Advance per-face flags and vaporized thickness at new timestep
    advanceTimeStep();

    // Total active plate height for Churchill-Chu natural convection
    const scalar plateHeight = max(sum(faceHeights_), SMALL);

    const scalarField alphaPrev = alpha_;

    forAll(h_, i)
    {
        const scalar p  = coolant_.coolantPressure()[i];
        const scalar T0 = coolant_.T0()[i];
        const scalar Tw = (*this)[i];

        const scalar Tsat         = IF97::Tsat97(p);
        const scalar TCHF         = Tsat + onbdTBessiron_;
        const scalar TLeidenfrost = Tsat + leidenfrostdTBessiron_;

        const bool postCrisis =
            vaporizedThicknessOld_[i] >= thicknessBessironRIA_;

        if (!postCrisis)
        {
            if (userTime < pulseStartTime_)
            {
                // Pre-pulse steady state: use single-phase NC and keep the
                // conduction clock at zero so it starts fresh at pulse onset.
                conductionTime_[i] = SMALL;
                h_[i] = singlePhase_.compute(T0, TCHF, p, 0.0, plateHeight);;
                alpha_[i] =
                    (1.0 - relax_)*alphaPrev[i]
                  + relax_*(alphaP[i] / (alphaP[i] + h_[i]));
            }
            else
            {
                // Transient conduction HTC — same formula for both sub-regimes.
                // conductionTime_ is zero before pulseStartTime_, so the clock
                // starts fresh at pulse onset.
                conductionTime_[i] =
                    max(conductionTimeOld_[i] + deltaT, SMALL);
                    
                const scalar k_liq       = IF97::tcond_Tp(T0, p);
                const scalar rho_liq     = IF97::rholiq_p(p);
                const scalar cp_liq      = IF97::cpmass_Tp(T0, p);
                const scalar alphaLiquid = k_liq / max(rho_liq*cp_liq, SMALL);

                const scalar hCond = 
                    k_liq
                  / sqrt
                    (
                        constant::mathematical::pi
                       *alphaLiquid
                       *conductionTime_[i]
                    );

                if (Tw < TCHF)
                {
                    // Pre-crisis: transient conduction in stagnant water
                    h_[i] = hCond;
                    alpha_[i] =
                        (1.0 - relax_)*alphaPrev[i]
                      + relax_*(alphaP[i] / (alphaP[i] + h_[i]));
                }
                else
                {
                    // Vaporizing: wall clamped at TCHF until crisis triggers.
                    // Use the heat flux at Tw = TCHF (not the previous-step heatFlux_,
                    // which may be elevated due to relaxation) so that the vaporization
                    // integral is consistent with the clamped wall temperature.
                    const scalar qAtTCHF =
                        alphaP[i] * (patchInternalField()()[i] - TCHF);

                    vaporizedThickness_[i] = computeVaporizedThickness
                    (
                        i, qAtTCHF, hCond, TCHF, T0, p, deltaT
                    );

                    if (vaporizedThickness_[i] >= thicknessBessironRIA_)
                    {
                        vaporizedThickness_[i] = thicknessBessironRIA_;
                        qCHF_[i] = qAtTCHF;
                    }

                    // Drive wall toward TCHF.
                    // h_equiv is the HTC that makes Twall = TCHF exactly:
                    //   alpha*T_int + (1-alpha)*T0 = TCHF
                    //   => h_equiv = alphaP*(T_int - TCHF)/(TCHF - T0)
                    h_[i] = max
                    (
                        alphaP[i]
                      * (patchInternalField()()[i] - TCHF)
                      / max(TCHF - T0, SMALL),
                        SMALL
                    );
                    alpha_[i] =
                        (1.0 - relax_)*alphaPrev[i]
                      + relax_*(alphaP[i] / (alphaP[i] + h_[i]));
                }
            }
        }
        else
        {
            // Initialise kBessiron_ once on the first post-crisis timestep.
            // Use the bare Sakurai curve (fbCorr=1) so k represents the
            // steady-state shape: as fbCorr decays to 1 the minimum of the
            // TB formula returns to exactly TLeidenfrost.
            if (kBessiron_[i] == 0.0)
            {
                const scalar qFilmLeid_0 =
                    filmBoiling_.compute
                    (
                        T0, TLeidenfrost, p, 0.0, 0.0, coolant_.DH()
                    )
                  * leidenfrostdTBessiron_;

                const scalar qFilmLeid_up_0 =
                    filmBoiling_.compute
                    (
                        T0, TLeidenfrost + 1.0, p, 0.0, 0.0, coolant_.DH()
                    )
                  * (leidenfrostdTBessiron_ + 1.0);

                kBessiron_[i] = computeTransitionDecay
                (
                    qCHF_[i],
                    qFilmLeid_0,
                    qFilmLeid_up_0 - qFilmLeid_0,
                    TLeidenfrost - TCHF
                );
            }

            if (Tw > Tsat)
            {
                if (Tw >= TCHF)
                {
                    // Film boiling correction: linearly decreases from initial
                    // value to 1 over filmBoilingRIACorrectionDuration_ s.
                    const scalar fbCorr = max
                    (
                        1.0,
                        filmBoilingRIACorrection_
                      - (filmBoilingRIACorrection_ - 1.0)
                      / filmBoilingRIACorrectionDuration_
                      * filmBoilingTimeOld_[i]
                    );

                    // qLeidenfrost stored for postprocessing only
                    qLeidenfrost_[i] =
                        filmBoiling_.compute
                        (
                            T0, TLeidenfrost, p, 0.0, 0.0, coolant_.DH()
                        )
                      * fbCorr * leidenfrostdTBessiron_;

                    // Unified post-crisis formula: Bessiron exponential with
                    // a live film boiling floor.  Above TLeidenfrost the
                    // exponential is negligible and the formula reduces to the
                    // corrected Sakurai film boiling heat flux.
                    const scalar qFilmTw =
                        filmBoiling_.compute
                        (
                            T0, Tw, p, 0.0, 0.0, coolant_.DH()
                        )
                      * fbCorr * (Tw - Tsat);

                    const scalar qTB =
                        qFilmTw
                      + (qCHF_[i] - qFilmTw)
                      * exp(-kBessiron_[i]*sqr(Tw - TCHF));

                    // Use (Tw - T0) as the driving delta-T so that h stays
                    // finite as Tw approaches Tsat (T0 < Tsat always).
                    h_[i] = qTB / max(Tw - T0, SMALL);

                    // fbCorr decays from the moment of crisis onward
                    filmBoilingTime_[i] = filmBoilingTimeOld_[i] + deltaT;
                }
                else
                {
                    // Tsat < Tw < TCHF: linear interpolation in q.
                    //
                    // The Bessiron formula gives q ≈ qCHF at Tw = TCHF, but
                    // still ≈ qCHF just above Tsat (exp factor ≈ 1 there),
                    // while NC gives a much smaller q at Tsat.  This
                    // discontinuity causes oscillation when Tw crosses Tsat.
                    // Bridging with a linear ramp makes q(Tw) continuous:
                    //   q(TCHF) = qCHF  (matches Bessiron from above)
                    //   q(Tsat) = q_NC(Tsat)  (matches NC from below)
                    const scalar hNC_atTsat =
                        singlePhase_.compute(T0, Tsat, p, 0.0, plateHeight);
                    const scalar qNC_atTsat =
                        hNC_atTsat * max(Tsat - T0, SMALL);

                    const scalar frac =
                        (Tw - Tsat) / max(TCHF - Tsat, SMALL);
                    const scalar qInterp =
                        qNC_atTsat + frac*(qCHF_[i] - qNC_atTsat);

                    h_[i] = qInterp / max(Tw - T0, SMALL);
                }

                alpha_[i] =
                    (1.0 - relax_)*alphaPrev[i]
                  + relax_*(alphaP[i] / (alphaP[i] + h_[i]));
            }
            else
            {
                // Natural convection after liquid contact
                h_[i] = singlePhase_.compute(T0, Tw, p, 0.0, plateHeight);

                alpha_[i] =
                    (1.0 - relax_)*alphaPrev[i]
                  + relax_*(alphaP[i] / (alphaP[i] + h_[i]));

                if (Tw <= T0 + resetDeltaTBessiron_)
                {
                    vaporizedThickness_[i] = 0.0;
                    conductionTime_[i] = 0.0;
                    filmBoilingTime_[i] = 0.0;
                    qCHF_[i] = 0.0;
                    qLeidenfrost_[i] = 0.0;
                    kBessiron_[i] = 0.0;
                }
            }
        }

        h_[i] = max(h_[i], SMALL);
    }

    fvPatchField<scalar>::updateCoeffs();
}


void coolantChannelRIAfvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    scalarField Tpred = alpha_*patchInternalField() + (1.0 - alpha_)*coolant_.T0();

    fvPatchField<scalar>::operator=(Tpred);

    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar>> coolantChannelRIAfvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar>> coolantChannelRIAfvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1.0 - alpha_)*coolant_.T0();
}


tmp<Field<scalar>> coolantChannelRIAfvPatchScalarField::gradientInternalCoeffs()
const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();
    return (alpha_ - 1.0)*deltaCoeffs;
}


tmp<Field<scalar>> coolantChannelRIAfvPatchScalarField::gradientBoundaryCoeffs()
const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();
    return (1.0 - alpha_)*deltaCoeffs*coolant_.T0();
}


void coolantChannelRIAfvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    coolant_.write(os);
    os.writeEntry("relax", relax_);
    os.writeEntry("correlations", correlations_);

    os.writeEntry("pulseStartTime",                    pulseStartTime_);
    os.writeEntry("filmBoilingRIACorrection",          filmBoilingRIACorrection_);
    os.writeEntry("filmBoilingRIACorrectionDuration",  filmBoilingRIACorrectionDuration_);
    os.writeEntry("leidenfrostdTBessiron",         leidenfrostdTBessiron_);
    os.writeEntry("onbdTBessiron",                 onbdTBessiron_);
    os.writeEntry("thicknessBessironRIA",          thicknessBessironRIA_);
    os.writeEntry("resetDeltaTBessiron",           resetDeltaTBessiron_);

    vaporizedThickness_.writeEntry("vaporizedThickness", os);
    qCHF_.writeEntry("qCHF", os);
    qLeidenfrost_.writeEntry("qLeidenfrost", os);
    kBessiron_.writeEntry("kBessiron", os);
    conductionTime_.writeEntry("conductionTime", os);
    filmBoilingTime_.writeEntry("filmBoilingTime", os);
    heatFlux_.writeEntry("heatFlux", os);
    h_.writeEntry("h", os);
    writeValueEntry(os);
}


makePatchTypeField(fvPatchScalarField, coolantChannelRIAfvPatchScalarField);

} // End namespace Foam

// ************************************************************************* //
