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

#include "coolantChannelfvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "globalOptions.H"
#include "wedgePolyPatch.H"
#include "IF97.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Private Functions  * * * * * * * * * * * * //

scalarField coolantChannelfvPatchScalarField::calculateFaceHeights() const
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
            const label pointIndex = pointLabels[pointI];
            const scalar z = localPoints[pointIndex] & pinDir;

            zMin = min(zMin, z);
            zMax = max(zMax, z);
        }

        faceHeights[faceI] = zMax - zMin;
    }

    return faceHeights;
}


void coolantChannelfvPatchScalarField::validateGeometry() const
{
    const polyBoundaryMesh& patches = patch().boundaryMesh().mesh().boundaryMesh();

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
            << "The coolant channel BC requires wedge patches so that the "
            << "channel surface represented by the wedge can be interpreted "
            << "consistently. No wedge patch was found."
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


void coolantChannelfvPatchScalarField::updateFaceHeightsIfNeeded()
{
    const scalar timeValue = this->db().time().value();

    if
    (
        timeValue > timeCheck_
     && db().foundObject<fvMesh>("referenceMesh")
    )
    {
        faceHeights_ = calculateFaceHeights();
    }
}


scalarField coolantChannelfvPatchScalarField::computeSolidSideConductance() const
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


void coolantChannelfvPatchScalarField::computeHeatFlux
(
    const scalarField& alphaP
)
{
    heatFlux_ = (patchInternalField() - *this) * alphaP;
}


scalar coolantChannelfvPatchScalarField::solvePreCHFTemperature
(
    const label faceI,
    const scalar qTarget
) const
{
    const scalar p = coolant_.coolantPressure()[faceI];
    const scalar T0 = coolant_.T0()[faceI];
    const scalar G = coolant_.massFlowRate()[faceI] / coolant_.A();
    const scalar hSP = singlePhase_.compute(T0, T0, p, G, coolant_.DH());

    // Bisect in [T0(or Tsat), Tmax] for T_w such that q_SB(T_w) = qTarget
    const scalar Tsat = IF97::Tsat97(p);
    scalar Tlo = max(T0, Tsat);
    scalar Thi = Tsat + 200.0;

    for (label iter = 0; iter < 20; ++iter)
    {
        const scalar Tmid = 0.5*(Tlo + Thi);

        const scalar hSB = (coolant_.Q()[faceI] <= 0.0)
            ? subcooledBoiling_.compute
              (
                  T0, Tmid, p, G, coolant_.DH(), qTarget, hSP
              )
            : saturatedBoiling_.compute
              (
                  T0, Tmid, p, coolant_.Q()[faceI], G, coolant_.DH(),
                  qTarget, hSP
              );

        const scalar qPred = hSB * (Tmid - T0);

        if (qPred < qTarget)
        {
            Tlo = Tmid;
        }
        else
        {
            Thi = Tmid;
        }

        if (mag(Thi - Tlo) < 0.01)
        {
            break;
        }
    }

    return 0.5*(Tlo + Thi);
}


scalar coolantChannelfvPatchScalarField::solveFilmBoilingTemperature
(
    const label faceI,
    const scalar qTarget
) const
{
    const scalar p = coolant_.coolantPressure()[faceI];
    const scalar T0 = coolant_.T0()[faceI];
    const scalar G = coolant_.massFlowRate()[faceI] / coolant_.A();
    const scalar Tsat = IF97::Tsat97(p);

    const scalar Tmax = Tsat + 2000.0;

    scalar Tlo = Tsat + 1e-3;
    scalar Thi = Tmax;

    for (label iter = 0; iter < 20; ++iter)
    {
        const scalar Tmid = 0.5*(Tlo + Thi);

        const scalar hFB =
            filmBoiling_.compute(T0, Tmid, p, coolant_.Q()[faceI], G, coolant_.DH());

        const scalar qPred = hFB * (Tmid - Tsat);

        if (qPred < qTarget)
        {
            Tlo = Tmid;
        }
        else
        {
            Thi = Tmid;
        }

        if (mag(Thi - Tlo) < 0.01)
        {
            break;
        }
    }

    return 0.5*(Tlo + Thi);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

coolantChannelfvPatchScalarField::coolantChannelfvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    kappaName_("k"),
    correlations_(),
    singlePhase_("Gnielinski"),
    onb_("Basu"),
    subcooledBoiling_("Rohsenow"),
    saturatedBoiling_("Rohsenow"),
    filmBoiling_("Frederking"),
    chf_("EPRI"),
    leidenfrost_("GroeneveldStewartCorrected"),
    transitionBoiling_("McDonoughMilichKing"),
    coolant_(p.size()),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    relax_(1.0),
    timeCheck_(0.0),
    heatFlux_(p.size(), 0.0),
    faceHeights_(p.size(), 0.0)
{}


coolantChannelfvPatchScalarField::coolantChannelfvPatchScalarField
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
    singlePhase_(correlations_.lookupOrDefault<word>("singlePhase", "Gnielinski")),
    onb_(correlations_.lookupOrDefault<word>("onb", "Basu")),
    subcooledBoiling_(correlations_.lookupOrDefault<word>("subcooledBoiling", "Rohsenow")),
    saturatedBoiling_(correlations_.lookupOrDefault<word>("saturatedBoiling", "Rohsenow")),
    filmBoiling_(correlations_.lookupOrDefault<word>("filmBoiling", "Frederking")),
    chf_(correlations_.lookupOrDefault<word>("chf", "EPRI")),
    leidenfrost_(correlations_.lookupOrDefault<word>("leidenfrost", "GroeneveldStewartCorrected")),
    transitionBoiling_(correlations_.lookupOrDefault<word>("transitionBoiling", "McDonoughMilichKing")),
    coolant_(dict, p.size(), chf_.model() == chfModel::Model::EPRI),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    timeCheck_(0.0),
    heatFlux_(p.size(), 0.0),
    faceHeights_(p.size(), 0.0)
{
    if (chf_.model() == chfModel::Model::GroeneveldTable)
    {
        chf_.initializeTable(dict, patch().boundaryMesh().mesh());
    }

    const scalar userTime =
        patch().boundaryMesh().mesh().time().timeToUserTime
        (
            this->db().time().value()
        );

    coolant_.updateBoundaryInputs(userTime);
    faceHeights_ = calculateFaceHeights();
    validateGeometry();
}


coolantChannelfvPatchScalarField::coolantChannelfvPatchScalarField
(
    const coolantChannelfvPatchScalarField& ptf,
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
    onb_(ptf.onb_),
    subcooledBoiling_(ptf.subcooledBoiling_),
    saturatedBoiling_(ptf.saturatedBoiling_),
    filmBoiling_(ptf.filmBoiling_),
    chf_(ptf.chf_),
    leidenfrost_(ptf.leidenfrost_),
    transitionBoiling_(ptf.transitionBoiling_),
    coolant_(ptf.coolant_),
    h_(ptf.h_),
    alpha_(ptf.alpha_),
    relax_(ptf.relax_),
    timeCheck_(ptf.timeCheck_),
    heatFlux_(ptf.heatFlux_),
    faceHeights_(ptf.faceHeights_)
{}


coolantChannelfvPatchScalarField::coolantChannelfvPatchScalarField
(
    const coolantChannelfvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(ptf, iF),
    kappaName_(ptf.kappaName_),
    correlations_(ptf.correlations_),
    singlePhase_(ptf.singlePhase_),
    onb_(ptf.onb_),
    subcooledBoiling_(ptf.subcooledBoiling_),
    saturatedBoiling_(ptf.saturatedBoiling_),
    filmBoiling_(ptf.filmBoiling_),
    chf_(ptf.chf_),
    leidenfrost_(ptf.leidenfrost_),
    transitionBoiling_(ptf.transitionBoiling_),
    coolant_(ptf.coolant_),
    h_(ptf.h_),
    alpha_(ptf.alpha_),
    relax_(ptf.relax_),
    timeCheck_(ptf.timeCheck_),
    heatFlux_(ptf.heatFlux_),
    faceHeights_(ptf.faceHeights_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coolantChannelfvPatchScalarField::updateCoeffs()
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

    const scalarField oldH = h_;

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

    coolant_.updateDerivedState(*this);

    // Pre-compute non-local patch quantities for CHF (no-op for non-EPRI models)
    chf_.precompute
    (
        zetas,
        faceHeights_,
        heatFlux_,
        coolant_.massFlowRate(),
        coolant_.coolantPressure(),
        coolant_.inletTemperature(),
        coolant_.A(),
        coolant_.flowDirection()
    );

    forAll(h_, i)
    {
        const scalar p  = coolant_.coolantPressure()[i];
        const scalar T0 = coolant_.T0()[i];
        const scalar Tw = (*this)[i];
        const scalar G  = coolant_.massFlowRate()[i] / coolant_.A();
        const scalar Q  = coolant_.Q()[i];

        const scalar hSP = singlePhase_.compute(T0, Tw, p, G, coolant_.DH());

        // Single-phase superheated vapour
        if (Q >= 0.999)
        {
            h_[i] = hSP;
            continue;
        }

        const scalar TONB = onb_.compute(T0, p, coolant_.coolantk()[i], hSP, i);

        // Single-phase liquid (below ONB temperature)
        if (Tw < TONB)
        {
            h_[i] = hSP;
            continue;
        }

        const scalar qCHF =
            chf_.compute(p, T0, Q, G, coolant_.DH(), heatFlux_[i], i);

        // Cheap predictor for the pre-CHF boiling branch at current wall temperature
        scalar hPreCHF = hSP;

        if (Q <= 0.0)
        {
            hPreCHF =
                subcooledBoiling_.compute
                (
                    T0, Tw, p, G, coolant_.DH(), heatFlux_[i], hSP
                );
        }
        else
        {
            hPreCHF =
                saturatedBoiling_.compute
                (
                    T0, Tw, p, Q, G, coolant_.DH(), heatFlux_[i], hSP
                );
        }

        // Cheap predictor for TCHF
        scalar TCHF = GREAT;

        if (mag(heatFlux_[i]) <= SMALL)
        {
            TCHF = GREAT;
        }
        else
        {
            TCHF = T0 + qCHF/max(hPreCHF, SMALL);
        }

        // Only solve near CHF. Can be tuned
        const scalar chfSolveBand = 0.10;

        if
        (
            qCHF < ROOTVGREAT
         && (
                mag(heatFlux_[i] - qCHF)/max(qCHF, SMALL) < chfSolveBand
             || mag(Tw - TCHF) < 5.0
            )
        )
        {
            TCHF = solvePreCHFTemperature(i, qCHF);
        }

        // Pre-CHF regime
        if (Tw < TCHF)
        {
            h_[i] = hPreCHF;
            h_[i] = max(h_[i], SMALL);
            continue;
        }

        const scalar hFB =
            filmBoiling_.compute(T0, Tw, p, Q, G, coolant_.DH());

        // T_Leidenfrost: direct for most correlations.
        // For Zuber (which gives q_Leid), solve q_FB(T_Leid) = q_Leid
        // self-consistently; only runs for post-CHF faces.
        const scalar TLeidenfrost = leidenfrost_.needsFilmBoilingTemperatureSolve()
            ? solveFilmBoilingTemperature(i, leidenfrost_.computeMinFlux(p))
            : leidenfrost_.compute(T0, p, Q, hFB);

        // Post-CHF: transition or film boiling
        if (Tw < TLeidenfrost)
        {
            h_[i] = transitionBoiling_.compute(p, heatFlux_[i], TCHF, qCHF, hFB);
        }
        else
        {
            h_[i] = hFB;
        }

        h_[i] = max(h_[i], SMALL);
    }

    if (sum(oldH) > 0.0)
    {
        h_ = (1.0 - relax_)*oldH + relax_*h_;
    }

    alpha_ = alphaP/(alphaP + h_);

    timeCheck_ = this->db().time().value();

    fvPatchField<scalar>::updateCoeffs();
}


void coolantChannelfvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<scalar>::operator=
    (
        alpha_*patchInternalField() + (1.0 - alpha_)*coolant_.T0()
    );

    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar>> coolantChannelfvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar>> coolantChannelfvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1.0 - alpha_)*coolant_.T0();
}


tmp<Field<scalar>> coolantChannelfvPatchScalarField::gradientInternalCoeffs() const
{
    return (alpha_ - 1.0)*patch().deltaCoeffs();
}


tmp<Field<scalar>> coolantChannelfvPatchScalarField::gradientBoundaryCoeffs() const
{
    return (1.0 - alpha_)*patch().deltaCoeffs()*coolant_.T0();
}


void coolantChannelfvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);

    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    coolant_.write(os);
    os.writeEntry("relax", relax_);
    os.writeEntry("correlations", correlations_);

    heatFlux_.writeEntry("heatFlux", os);
    h_.writeEntry("h", os);
    writeValueEntry(os);
}


makePatchTypeField(fvPatchScalarField, coolantChannelfvPatchScalarField);

} // End namespace Foam

// ************************************************************************* //
