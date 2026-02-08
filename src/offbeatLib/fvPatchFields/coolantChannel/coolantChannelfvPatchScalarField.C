/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "coolantChannelfvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "scalarFieldFieldINew.H"
#include "PtrListScalarFieldFieldINew.H"
#include "IF97.H"
#include "globalOptions.H"
#include "wedgePolyPatch.H"
#include "Pstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Private functions * * * * * * * * * * * * * * * * //

Foam::scalarField coolantChannelfvPatchScalarField::calculateFaceHeights()
{
    const fvPatch& patch = this->patch();
    const pointField& localPoints = patch.patch().localPoints();
    const globalOptions& globalOpt
    (
        patch.patch().boundaryMesh().mesh().lookupObject<globalOptions>
        ("globalOptions")
    );
    const vector pinDirection(globalOpt.pinDirection());

    scalarField faceHeights(this->size(), 0.0);

    forAll(patch.patch().localFaces(), faceI)
    {
        const labelList& pointLabelList = patch.patch().localFaces()[faceI];
        scalar maxHeight = -GREAT;
        scalar minHeight = GREAT;

        forAll(pointLabelList, pointI)
        {
            const label& pointIndex = pointLabelList[pointI];
            
            maxHeight = max
            (
                maxHeight, 
                localPoints[pointIndex] & pinDirection
            );
            
            minHeight = min
            (
                minHeight, 
                localPoints[pointIndex] & pinDirection
            );
        }

        faceHeights[faceI] = maxHeight - minHeight; 
    }

    return faceHeights;
}

Foam::dictionary coolantChannelfvPatchScalarField::initCorrelationsDict
(
    const dictionary& dict
)
{
    dictionary d = dict.subOrEmptyDict("correlations");

    // Set default correlations
    if(!d.found("ONB")) d.add("ONB", "Basu");
    if(!d.found("CHF")) d.add("CHF", "EPRI");
    if(!d.found("Leidenfrost")) d.add("Leidenfrost", "GroeneveldStewartCorrected");
    if(!d.found("SinglePhase")) d.add("SinglePhase", "Gnielinski");
    if(!d.found("SubcooledBoiling")) d.add("SubcooledBoiling", "Rohsenow");
    if(!d.found("SaturatedBoiling")) d.add("SaturatedBoiling", "Rohsenow");
    if(!d.found("TransitionBoiling")) d.add("TransitionBoiling", "McDonoughMilichKing");
    if(!d.found("FilmBoiling")) d.add("FilmBoiling", "Frederking");

    return d;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coolantChannelfvPatchScalarField::
coolantChannelfvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(p, iF),
    kappaName_("k"),
    correlations_(),
    correlationRegimeList_(p.size(), "empty"),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    A_(-1.0),
    flowDirection_(1.0),
    relaxHTC_(0.15),
    filmBoilingRIACorrection_(1.0),
    transitionBoilingRIACorrection_(-0.00006),
    leidenfrostdTBessiron_(50.0),
    thicknessBessironRIA_(0.000160),
    timeCheck_(0.0),
    DH_(-1.0),
    qCHF_(p.size(), 0.0),
    vaporizedThickness_(p.size(), 0.0),
    heatFlux_(p.size(), 0.0),
    Q_(p.size(), 0.0),
    X_(p.size(), 0.0),
    Re_(p.size(), 0.0),
    Pr_(p.size(), 0.0),
    Twall_(p.size(), 0.0),
    Tfilm_(p.size(), 0.0),
    coolantk_(p.size(), 0.0),
    coolantEnthalpy_(p.size(), 0.0),
    coolantPressure_(p.size(), 0.0),
    coolantPressureList_(),
    massFlowRate_(p.size(), -1.0),
    massFlowRateList_(),
    inletTemperature_(-1.0),
    inletTemperatureList_(),
    enthalpyModel_(),
    zValues_(),
    timeValues_(),
    T0ProfileData_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    dict_(),
    faceHeights_(p.size(), 0.0),
    PsiToPa_(6894.757),
    BtuToJ_(1055.06),
    FeetToMeters_(0.3048),
    LbToKg_(0.45359237),
    HrToSec_(3600.0),
    g_(9.80665)
{}


coolantChannelfvPatchScalarField::
coolantChannelfvPatchScalarField
 (
     const fvPatch& p,
     const DimensionedField<scalar, volMesh>& iF,
     const dictionary& dict,
     const bool valueRequired
 )
 :
    fvPatchField<scalar>(p, iF, dict, valueRequired),
    kappaName_
    (
        dict.lookupOrDefault<word>("kappa", "k")
    ),
    correlations_(initCorrelationsDict(dict)),
    correlationRegimeList_(this->size(), "empty"),
    T0_(p.size(), 0.0),
    h_(p.size(), 0.0),
    alpha_(p.size(), 0.0),
    A_(readScalar(dict.lookup("flowArea"))),
    flowDirection_(dict.lookupOrDefault("flowDirection", 1.0)),
    relaxHTC_
    (
        dict.lookupOrDefault("relaxHTC", 0.15)
    ),
    filmBoilingRIACorrection_
    (
        dict.lookupOrDefault("filmBoilingRIACorrection", 1.0)
    ),
    transitionBoilingRIACorrection_
    (
        dict.lookupOrDefault("transitionBoilingRIACorrection", -0.00006)
    ),
    leidenfrostdTBessiron_
    (
        dict.lookupOrDefault("leidenfrostdTBessiron", 50.0)
    ),
    thicknessBessironRIA_
    (
        dict.lookupOrDefault("thicknessBessironRIA", 0.000160)
    ),
    timeCheck_(0.0),
    DH_(readScalar(dict.lookup("hydraulicDiameter"))),
    qCHF_(p.size(), 0.0),
    vaporizedThickness_(p.size(), 0.0),
    heatFlux_(p.size(), 0.0),
    Q_(p.size(), 0.0),
    X_(p.size(), 0.0),
    Re_(p.size(), 0.0),
    Pr_(p.size(), 0.0),
    Twall_(p.size(), 0.0),
    Tfilm_(p.size(), 0.0),
    coolantk_(p.size(), 0.0),
    coolantEnthalpy_(p.size(), 0.0),
    coolantPressure_(p.size(), 0),
    coolantPressureList_(),
    massFlowRate_(p.size(), -1.0),
    massFlowRateList_(),
    inletTemperature_(-1.0),
    inletTemperatureList_(),
#ifdef OPENFOAMFOUNDATION
    enthalpyModel_(dict.lookup<bool>("enthalpyModel")),
#elif OPENFOAMESI
    enthalpyModel_(dict.get<bool>("enthalpyModel")),
#endif
    PsiToPa_(6894.757),
    BtuToJ_(1055.06),
    FeetToMeters_(0.3048),
    LbToKg_(0.45359237),
    HrToSec_(3600.0),
    g_(9.80665)
{
    if
    (
        dict.found("qCHF") 
        and 
        correlations_.lookup<word>("CHF") == "Bessiron"
    )
    {
        qCHF_ = scalarField("qCHF", dict, p.size());
    }

    if
    (
        dict.found("vaporizedThickness") 
        and 
        correlations_.lookup<word>("CHF") == "Bessiron"
    )
    {
        vaporizedThickness_ = scalarField("vaporizedThickness", dict, p.size());
    }
          
    if (correlations_.lookup<word>("CHF") == "lookuptableGroeneveld")
    {
        // Instantiate IO dictionary for CHF data
        IOdictionary chfDict
        (
            IOobject
            (
                "chfTableDict",
                patch().boundaryMesh().mesh().time().constant(),
                patch().boundaryMesh().mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        chfData_ = PtrList<FieldField<Field, scalar>>
            (   
                PtrList<FieldField<Field, scalar>>(chfDict.lookup("data"), 
                PtrListScalarFieldFieldINew())
            );
            
        chfPitchToDiameter_ = scalar
        (
            chfDict.lookupOrDefault<scalar>("PitchToDiameter", 0.0)
        );
        chfPressureValues_ = scalarField(chfDict.lookup("pressureValues"));
        chfMassFluxValues_ = scalarField(chfDict.lookup("massFluxValues"));
        chfQualityValues_ = scalarField(chfDict.lookup("qualityValues"));
        chfPressureMethod_ = interpolateTableBase::interpolationMethod
            (
                interpolateTableBase::interpolationMethodNames_
                [
                    dict.lookupOrDefault<word>("pressureInterpolationMethod", "linear")
                ]
            );
        chfMassFluxMethod_ = interpolateTableBase::interpolationMethod
            (
                interpolateTableBase::interpolationMethodNames_
                [
                    dict.lookupOrDefault<word>("massFluxInterpolationMethod", "linear")
                ]
            );
        chfQualityMethod_ = interpolateTableBase::interpolationMethod
            (interpolateTableBase::interpolationMethodNames_
                [
                    dict.lookupOrDefault<word>("qualityInterpolationMethod", "linear")
                ]
            );
        
#ifdef OPENFOAMFOUNDATION            
        chfPressureTable_.set
        ( 
            new scalarFieldFieldInterpolateTable
            (chfPressureValues_, chfData_, chfPressureMethod_)
        );
#elif OPENFOAMESI
        chfPressureTable_.reset
        (
            new scalarFieldFieldInterpolateTable
            (chfPressureValues_, chfData_, chfPressureMethod_)  
        );
#endif 
    }

    if(enthalpyModel_ == false)
    {
        T0ProfileData_ = FieldField<Field, scalar>
            (
                PtrList<scalarField>(dict.subDict("axialProfileDict").lookup("T0Data"), 
                    scalarFieldFieldINew())
            );

        zMethod_ = interpolateTableBase::interpolationMethod
            (
                interpolateTableBase::interpolationMethodNames_
                [
                    dict.subDict("axialProfileDict").lookupOrDefault<word>
                    ("axialInterpolationMethod", "linear")
                ]
            );
        tMethod_ = interpolateTableBase::interpolationMethod
            (
                interpolateTableBase::interpolationMethodNames_
                [
                    dict.subDict("axialProfileDict").lookupOrDefault<word>
                    ("timeInterpolationMethod", "linear")
                ]
            );
        zValues_ = scalarField(dict.subDict("axialProfileDict").lookup("axialLocations"));
        timeValues_ = scalarField(dict.subDict("axialProfileDict").lookup("timePoints"));
        dict_ = dictionary(dict.subDict("axialProfileDict"));

#ifdef OPENFOAMFOUNDATION            
        T0Table_.set
        ( 
            new scalarFieldInterpolateTable
            (timeValues_, T0ProfileData_, tMethod_)
        );
#elif OPENFOAMESI
        T0Table_.reset
        (
            new scalarFieldInterpolateTable
            (timeValues_, T0ProfileData_, tMethod_)  
        );
#endif 
    }

    //- Read fluid pressure (either list or fixed pressure)
    if(dict.found("coolantPressureList"))
    {
#ifdef OPENFOAMFOUNDATION            
        coolantPressureList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "coolantPressureList", dict.subDict("coolantPressureList")
            )
        );
#elif OPENFOAMESI            
        coolantPressureList_.reset
        ( 
            new scalarTable
            (
                "coolantPressureList", dict.subDict("coolantPressureList")
            )
        );
#endif        
        
        coolantPressure_ = 
        coolantPressureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }
    else if (dict.found("coolantPressure"))
    {
        coolantPressure_ = scalarField("coolantPressure", dict, p.size());
    }
    else
    {
        FatalErrorInFunction() 
            << "Coolant channel BC for " << patch().name() << " requires:" << nl
            << "- \"coolantPressureList\" or" << nl
            << "- \"coolantPressure\"" << abort(FatalError) << endl;
    }

    //- Read fluid mass flow rate (either list or fixed)
    if(dict.found("massFlowRateList"))
    {
#ifdef OPENFOAMFOUNDATION            
        massFlowRateList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "massFlowRateList", dict.subDict("massFlowRateList")
            )
        );
#elif OPENFOAMESI            
        massFlowRateList_.reset
        ( 
            new scalarTable
            (
                "massFlowRateList", dict.subDict("massFlowRateList")
            )
        );
#endif        
        
        massFlowRate_ = 
        massFlowRateList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }
    else if (dict.found("massFlowRate"))
    {
        massFlowRate_ = scalarField("massFlowRate", dict, p.size());
    }
    else
    {
        FatalErrorInFunction() 
            << "Coolant channel BC for " << patch().name() << " requires:" << nl
            << "- \"massFlowRateList\" or" << nl
            << "- \"massFlowRate\"" << abort(FatalError) << endl;
    }

    if (enthalpyModel_ == true || correlations_.lookup<word>("CHF") == "EPRI")
    {
        //- Read fluid inlet temperature (either list or fixed)
        if(dict.found("inletTemperatureList"))
        {
#ifdef OPENFOAMFOUNDATION            
            inletTemperatureList_.set
            ( 
                new Function1s::Table<scalar>
                (
                    "inletTemperatureList", dict.subDict("inletTemperatureList")
                )
            );
#elif OPENFOAMESI            
            inletTemperatureList_.reset
            ( 
                new scalarTable
                (
                    "inletTemperatureList", dict.subDict("inletTemperatureList")
                )
            );
#endif        
            
            inletTemperature_ = 
            inletTemperatureList_->value
            (
                patch().boundaryMesh().mesh().time().timeToUserTime
                (
                    this->db().time().value()
                )
            );
        }
        else if (dict.found("inletTemperature"))
        {
            inletTemperature_ = scalar(readScalar(dict.lookup("inletTemperature")));
        }
        else
        {
            FatalErrorInFunction() 
                << "Coolant channel BC for " << patch().name() << " requires:" << nl
                << "- \"inletTemperatureList\" or" << nl
                << "- \"inletTemperature\"" << abort(FatalError) << endl;
        }
    }

    faceHeights_ = scalarField(calculateFaceHeights());

    const polyBoundaryMesh& patches = patch().boundaryMesh().mesh().boundaryMesh();
    bool checkWedge(false);
    forAll(patch().boundaryMesh(), i)
    {
        if
        (
            patches[i].size() && isType<wedgePolyPatch>(patches[i])
        )
        {
            checkWedge = true;
        }
    }
    if (checkWedge == false)
    {
        FatalErrorIn("coolantChannelfvPatchScalarField")
            << "The model does not contain wedge patches and automatic " << nl
            << "calculation of model angularFraction is not possible." << nl
            << "Please provide angularFraction keyword in the globalOptions dict"
            << exit(FatalError);
    }
}

coolantChannelfvPatchScalarField::
coolantChannelfvPatchScalarField
(
    const coolantChannelfvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
#ifdef OPENFOAMFOUNDATION
    fvPatchField<scalar>(ptf, p, iF, mapper, mappingRequired),
#elif OPENFOAMESI
    fvPatchField<scalar>(ptf, p, iF, mapper),
#endif
    kappaName_(ptf.kappaName_),
    correlations_(ptf.correlations_),
    correlationRegimeList_(ptf.correlationRegimeList_),
    T0_(ptf.T0_),
    h_(ptf.h_),
    alpha_(ptf.alpha_),
    A_(ptf.A_),
    flowDirection_(ptf.flowDirection_),
    relaxHTC_(ptf.relaxHTC_),
    timeCheck_(ptf.timeCheck_),
    DH_(ptf.DH_),
    qCHF_(ptf.qCHF_),
    vaporizedThickness_(ptf.vaporizedThickness_),
    heatFlux_(ptf.heatFlux_),
    Q_(ptf.Q_),
    X_(ptf.X_),
    Re_(ptf.Re_),
    Pr_(ptf.Pr_),
    Twall_(ptf.Twall_),
    Tfilm_(ptf.Tfilm_),
    coolantk_(ptf.coolantk_),
    coolantEnthalpy_(ptf.coolantEnthalpy_),
    coolantPressure_(ptf.coolantPressure_),
    coolantPressureList_(ptf.coolantPressureList_),
    massFlowRate_(ptf.massFlowRate_),
    massFlowRateList_(ptf.massFlowRateList_),
    inletTemperature_(ptf.inletTemperature_),
    inletTemperatureList_(ptf.inletTemperatureList_),
    enthalpyModel_(ptf.enthalpyModel_),
    zValues_(ptf.zValues_),
    timeValues_(ptf.timeValues_),
    T0ProfileData_(ptf.T0ProfileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    dict_(),
    PsiToPa_(ptf.PsiToPa_),
    BtuToJ_(ptf.BtuToJ_),
    FeetToMeters_(ptf.FeetToMeters_),
    LbToKg_(ptf.LbToKg_),
    HrToSec_(ptf.HrToSec_),
    g_(ptf.g_)
{}


coolantChannelfvPatchScalarField::
coolantChannelfvPatchScalarField
(
    const coolantChannelfvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fvPatchField<scalar>(tppsf, iF),
    kappaName_(tppsf.kappaName_),
    correlations_(tppsf.correlations_),
    correlationRegimeList_(tppsf.correlationRegimeList_),
    T0_(tppsf.T0_),
    h_(tppsf.h_),
    alpha_(tppsf.alpha_),
    A_(tppsf.A_),
    flowDirection_(tppsf.flowDirection_),
    relaxHTC_(tppsf.relaxHTC_),
    timeCheck_(tppsf.timeCheck_),
    DH_(tppsf.DH_),
    qCHF_(tppsf.Q_),
    vaporizedThickness_(tppsf.vaporizedThickness_),
    heatFlux_(tppsf.heatFlux_),
    Q_(tppsf.Q_),
    X_(tppsf.X_),
    Re_(tppsf.Re_),
    Pr_(tppsf.Pr_),
    Twall_(tppsf.Twall_),
    Tfilm_(tppsf.Tfilm_),
    coolantk_(tppsf.coolantk_),
    coolantEnthalpy_(tppsf.coolantEnthalpy_),
    coolantPressure_(tppsf.coolantPressure_),
    coolantPressureList_(),
    massFlowRate_(tppsf.massFlowRate_),
    massFlowRateList_(),
    inletTemperature_(tppsf.inletTemperature_),
    inletTemperatureList_(),
    enthalpyModel_(tppsf.enthalpyModel_),
    zValues_(tppsf.zValues_),
    timeValues_(tppsf.timeValues_),
    T0ProfileData_(tppsf.T0ProfileData_),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    dict_(),
    PsiToPa_(tppsf.PsiToPa_),
    BtuToJ_(tppsf.BtuToJ_),
    FeetToMeters_(tppsf.FeetToMeters_),
    LbToKg_(tppsf.LbToKg_),
    HrToSec_(tppsf.HrToSec_),
    g_(tppsf.g_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * *s * * * * * * * //

void coolantChannelfvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Update face heights in every timestep for DD solvers
    if (this->patch().boundaryMesh().mesh().time().value() > timeCheck_ &&
        db().foundObject<fvMesh>("referenceMesh"))
    {
        faceHeights_ = scalarField(calculateFaceHeights());
    }

    // Update fluid pressure
    if (coolantPressureList_.valid())
    {
        coolantPressure_ = 
        coolantPressureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    // Update mass flow rate
    if (massFlowRateList_.valid())
    {
        massFlowRate_ = 
        massFlowRateList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    // Update inlet temperature
    if (inletTemperatureList_.valid())
    {
        inletTemperature_ = 
        inletTemperatureList_->value
        (
            patch().boundaryMesh().mesh().time().timeToUserTime
            (
                this->db().time().value()
            )
        );
    }

    // Compute face radii and axial elevation of the face centers
    const globalOptions& globalOpt
    (
        patch().boundaryMesh().mesh().lookupObject<globalOptions>
        ("globalOptions")
    );
    const vector pinDirection(globalOpt.pinDirection());
    scalarField zetas = patch().Cf() & pinDirection;
    scalarField radii = mag(patch().Cf() - zetas * pinDirection);
    scalar minZeta = min(zetas);
    reduce(minZeta, minOp<scalar>());

    // Updating the wall temperature
    Twall_ = *this;

    // Initializing fuel conductivity and wall heat flux
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const fvPatch& patch = this->patch();
    const scalarField& deltaCoeffs = patch.deltaCoeffs();
    scalarField alphaN(this->size(), 0.0);
    scalarField alphaP(this->size(), 0.0);

    if (mesh.foundObject<volTensorField>(kappaName_))
    {
        const fvPatchField<tensor>& k = 
            patch.lookupPatchField<volTensorField, tensor>
            (
                kappaName_
            );
        
        vectorField n = (patch.Sf() / patch.magSf());
        alphaP = ((n & k) & n) * deltaCoeffs;
        heatFlux_ = (patchInternalField() - Twall_) * alphaP;
    }
    else
    {
        const fvPatchField<scalar>& k = 
            patch.lookupPatchField<volScalarField, scalar>
            (
                kappaName_
            );

        alphaP = k * deltaCoeffs;
        heatFlux_ = (patchInternalField() - Twall_) * alphaP;
    }

    scalarField oldH(this->size(), 0.0);
    scalarField beta(this->size(), 0.0);

    // Preparing enthalpy increments for the enthalpy model
    scalarField combinedZetas;
    scalarField combinedEnthalpyIncrements;
    if (enthalpyModel_ == true)
    {
        if (Pstream::parRun())
        {
            // Gathering and combining enthalpy increments for parallel cases
            List<scalarField> faceZetas(Pstream::nProcs(), zetas);
            List<scalarField> faceEnthalpyIncrements(Pstream::nProcs(), 2.0 * radii * 
                constant::mathematical::pi * heatFlux_ * faceHeights_ / massFlowRate_);
            Pstream::gatherList(faceZetas);
            Pstream::gatherList(faceEnthalpyIncrements);
            Pstream::scatterList(faceZetas);
            Pstream::scatterList(faceEnthalpyIncrements);
            forAll(faceZetas, i)
            {
                combinedZetas.append(faceZetas[i]);
                combinedEnthalpyIncrements.append(faceEnthalpyIncrements[i]);
            }
        }
        else
        {
            combinedZetas = zetas;
            combinedEnthalpyIncrements = 2.0 * radii * 
                constant::mathematical::pi * heatFlux_ * faceHeights_ / massFlowRate_;
        }
    }

    // Calculating new coolant TH profiles
    forAll(zetas, i)
    {
        if (enthalpyModel_ == true)
        {
            coolantEnthalpy_[i] = IF97::hmass_Tp(inletTemperature_, coolantPressure_[i]);

            forAll(combinedZetas, j)
            {
                if (combinedZetas[j] * sign(flowDirection_) <= zetas[i] * sign(flowDirection_))
                {
                    coolantEnthalpy_[i] = coolantEnthalpy_[i] + combinedEnthalpyIncrements[j];
                }
            }

            T0_[i] = IF97::T_phmass(coolantPressure_[i], coolantEnthalpy_[i]);
        }
        else
        {
            scalarField T0_zData
            (
                T0Table_()
                (
                    this->patch().boundaryMesh().mesh().time().timeToUserTime
                    (
                        this->db().time().value()
                    )
                )
            );
            scalarInterpolateTable T0_zTable(zValues_, T0_zData, zMethod_, 
                interpolateTableBase::outOfBoundsMethod::FIXED);
            T0_[i] = T0_zTable(zetas[i]);
            coolantEnthalpy_[i] = IF97::hmass_Tp(T0_[i], coolantPressure_[i]);
        }
        
        Q_[i] = (coolantEnthalpy_[i] - IF97::hliq_p(coolantPressure_[i])) / 
            (IF97::hvap_p(coolantPressure_[i]) - IF97::hliq_p(coolantPressure_[i]));
        Tfilm_[i] = (T0_[i] + Twall_[i]) * 0.5;
        
        if (Q_[i] <= 0.0)
        {
            X_[i] = 0.0;
        }
        else if (Q_[i] >= 1.0)
        {
            X_[i] = 1.0;
        }
        else
        {
            // Calculating void fraction assuming the slip ratio to be unity (HEM model).
            X_[i] = 1.0 / (1.0 + (1.0 - Q_[i]) / Q_[i] * 
                IF97::rhovap_p(coolantPressure_[i]) / IF97::rholiq_p(coolantPressure_[i]));
        }

        Re_[i] = massFlowRate_[i] * DH_ / IF97::visc_Tp(Tfilm_[i], coolantPressure_[i]) / A_;
        Pr_[i] = IF97::prandtl_Tp(Tfilm_[i], coolantPressure_[i]);
        coolantk_[i] = IF97::tcond_Tp(Tfilm_[i], coolantPressure_[i]);
        oldH = h_;

        // Volumetric thermal expansion coefficient: Temporary table added to IF97.H
        // Valid for liquid water at atmospheric pressure (278.15 K < T < 373.15 K, p = 101325 Pa)
        // Temperature values beyond temperature limits are clamped.
        beta[i] = IF97::volumetricThermalExpansion_T(Tfilm_[i]);
    }

    // Empty HTC vectors for the different boiling regimes
    scalarField hSP(this->size(), 0.0);
    scalarField hSCB(this->size(), 0.0);
    scalarField hSB(this->size(), 0.0);
    scalarField hFB(this->size(), 0.0);
    scalarField hTB(this->size(), 0.0);

    // Calculating new heat transfer coefficients
    forAll(zetas, i)
    {
        // Single phase regime HTC
        // Analytical value (laminar regime, Re < 2300)
        if (correlations_.lookup<word>("SinglePhase") == "analyticalLaminar")
        {
            hSP[i] = coolantk_[i] / DH_ * 4.36;
        }
        // Gnielinski correlation with Filonenko friction factor 
        // (1976, transition and turbulent regimes; 2300 < Re < 5e6, 0.5 < Pr < 2000)
        else if (correlations_.lookup<word>("SinglePhase") == "Gnielinski")
        {
            scalar filonenko(pow(1.58 * log(Re_[i]) - 3.28, -2.0));
            hSP[i] = coolantk_[i] / DH_ * ((0.5 * filonenko) * (Re_[i] - 1000.0) * Pr_[i]) / 
                (1.0 + 12.7 * sqrt(0.5 * filonenko) * (pow(Pr_[i], 2.0 / 3.0) - 1.0));
        }
        // Natural convection HTC correlation for tube geometry (TRACE)
        else if (correlations_.lookup<word>("SinglePhase") == "TRACETubeNC")
        {
            scalar hNCturb_(0.1 * coolantk_[i] / DH_ * pow(Pr_[i], 1.0 / 3.0) * 
                cbrt((g_ * beta[i]* (Twall_[i] - T0_[i]) * pow(DH_, 3.0) * 
                pow(IF97::rhomass_Tp(Tfilm_[i], coolantPressure_[i]), 2.0)) / 
                pow(IF97::visc_Tp(Tfilm_[i], coolantPressure_[i]), 2.0)));
            scalar hNClam_(0.59 * coolantk_[i] / DH_ * pow(Pr_[i], 0.25) * 
                pow((g_ * beta[i] * mag(Twall_[i] - T0_[i]) * pow(DH_, 3.0) * 
                pow(IF97::rhomass_Tp(Tfilm_[i], coolantPressure_[i]), 2.0)) / 
                pow(IF97::visc_Tp(Tfilm_[i], coolantPressure_[i]), 2.0), 0.25));
            hSP[i] = max(hNCturb_, hNClam_);
        }
        // Natural convection HTC correlation for rod bundle (TRACE)
        else if (correlations_.lookup<word>("SinglePhase") == "TRACERodBundleNC")
        {
            hSP[i] = 0.7 * coolantk_[i] / DH_ * pow(Pr_[i], 0.25) * 
                pow((g_ * beta[i] * mag(Twall_[i] - T0_[i]) * pow(DH_, 3.0) * 
                pow(IF97::rhomass_Tp(Tfilm_[i], coolantPressure_[i]), 2.0)) / 
                pow(IF97::visc_Tp(Tfilm_[i], coolantPressure_[i]), 2.0), 0.25);
        }
        // Dittus-Boelter correlation 
        // (1930, legacy; only turbulent regime; 1e4 < Re < 1.24e5, 0.6 < Pr < 160)
        else if (correlations_.lookup<word>("SinglePhase") == "DittusBoelter")
        {
            hSP[i] = coolantk_[i] / DH_ * 0.023 * pow(Re_[i], 0.8) * pow(Pr_[i], 0.4);
        }
        else
        {
            FatalErrorInFunction()  
                << "Incorrectly specified single phase HTC correlation name."  
                << abort(FatalError);
        }
        
        // Subcooled boiling regime HTC
        // Chen's correlation (1963)
        if (correlations_.lookup<word>("SubcooledBoiling") == "Chen")
        {
            scalar F(1.0);

            scalar hc(0.023 * pow((massFlowRate_[i] / A_ * DH_) / 
                IF97::viscliq_p(coolantPressure_[i]), 0.8) * 
                pow(IF97::prandtlliq_p(coolantPressure_[i]), 0.4) * 
                IF97::tcondliq_p(coolantPressure_[i]) / DH_);

            scalar S(1.0 / (1.0 + 2.53 * 1e-6 * pow(massFlowRate_[i] / A_ * DH_ / 
                IF97::visc_Tp(Tfilm_[i], coolantPressure_[i]) * pow(F, 1.25), 1.17)));

            scalar hNB(0.00122 * pow(Twall_[i] - IF97::Tsat97(coolantPressure_[i]), 0.24) * 
                pow(IF97::psat97(Twall_[i]) - 
                IF97::psat97(IF97::Tsat97(coolantPressure_[i])), 0.75) * 
                (pow(IF97::tcondliq_p(coolantPressure_[i]), 0.79) * 
                pow(IF97::cpliq_p(coolantPressure_[i]), 0.45) * 
                pow(IF97::rholiq_p(coolantPressure_[i]), 0.49)) / 
                (pow(IF97::sigma97(IF97::Tsat97(coolantPressure_[i])), 0.5) * 
                pow(IF97::viscliq_p(coolantPressure_[i]), 0.29) * 
                pow(IF97::latentHeatVaporization_p(coolantPressure_[i]), 0.24) *
                pow(IF97::rhovap_p(coolantPressure_[i]), 0.24)));

            hSCB[i] = F * hc + S * hNB;
        }
        // Jens-Lottes correlation (1951)
        else if (correlations_.lookup<word>("SubcooledBoiling") == "JensLottes")
        {
            hSCB[i] = heatFlux_[i] / (25.0 * pow(mag(heatFlux_[i]) / 1e6, 0.25) / 
                exp(coolantPressure_[i] / 6200000.0) + 
                IF97::Tsat97(coolantPressure_[i]) - T0_[i]);
        }
        // Thom's correlation (1965)
        else if (correlations_.lookup<word>("SubcooledBoiling") == "Thom")
        {
            hSCB[i] = heatFlux_[i] / (22.65 * pow(mag(heatFlux_[i]) / 1e6, 0.5) / 
                exp(coolantPressure_[i] / 8700000.0) + 
                IF97::Tsat97(coolantPressure_[i]) - T0_[i]);
        }
        // Rohsenow correlation (1952)
        else if (correlations_.lookup<word>("SubcooledBoiling") == "Rohsenow")
        {
            hSCB[i] = heatFlux_[i] / (cbrt(heatFlux_[i] / 
                (IF97::viscliq_p(coolantPressure_[i]) * 
                IF97::latentHeatVaporization_p(coolantPressure_[i]) * 
                pow(g_ * (IF97::rholiq_p(coolantPressure_[i]) - 
                IF97::rhovap_p(coolantPressure_[i])) / 
                IF97::sigma97(IF97::Tsat97(coolantPressure_[i])), 0.5))) *
                (0.013 * IF97::latentHeatVaporization_p(coolantPressure_[i]) * 
                pow(IF97::prandtlliq_p(coolantPressure_[i]), 1.0)) / 
                IF97::cpliq_p(coolantPressure_[i]) +
                IF97::Tsat97(coolantPressure_[i]) - T0_[i]);
        }
        else
        {
            FatalErrorInFunction()  
                << "Incorrectly specified subcooled boiling HTC correlation name."  
                << abort(FatalError);
        }

        // Saturated boiling regime HTC
        // Chen's correlation (1963)
        if (correlations_.lookup<word>("SaturatedBoiling") == "Chen")
        {
            scalar tempQ(0.0);
            if (Q_[i] < 0.0)
            {
                tempQ = 0.0;
            }
            else
            {
                tempQ = Q_[i];
            }
            scalar F(0.0);
            scalar recXtt(pow(tempQ / (1 - tempQ), 0.9) * 
                pow(IF97::rholiq_p(coolantPressure_[i]) / 
                IF97::rhovap_p(coolantPressure_[i]), 0.5) * 
                pow(IF97::viscvap_p(coolantPressure_[i]) / 
                IF97::viscliq_p(coolantPressure_[i]), 0.1));

            scalar hc(0.023 * pow((massFlowRate_[i] / A_ * DH_ * (1.0 - tempQ)) / 
                IF97::viscliq_p(coolantPressure_[i]), 0.8) * 
                pow(IF97::prandtlliq_p(coolantPressure_[i]), 0.4) * 
                IF97::tcondliq_p(coolantPressure_[i]) / DH_);

            if (recXtt < 0.1)
            {
                F = 1.0;
            }
            else
            {
                F = 2.35 * pow(0.213 + recXtt, 0.736);
            }

            scalar S(1.0 / (1.0 + 2.53 * 1e-6 * pow(massFlowRate_[i] / A_ * DH_ / 
                IF97::visc_Tp(Tfilm_[i], coolantPressure_[i]) * (1.0 - Q_[i]) * pow(F, 1.25), 1.17)));

            scalar hNB(0.00122 * pow(mag(Twall_[i] - IF97::Tsat97(coolantPressure_[i])), 0.24) * 
                pow(mag(IF97::psat97(Twall_[i]) - 
                IF97::psat97(IF97::Tsat97(coolantPressure_[i]))), 0.75) * 
                (pow(IF97::tcondliq_p(coolantPressure_[i]), 0.79) * 
                pow(IF97::cpliq_p(coolantPressure_[i]), 0.45) * 
                pow(IF97::rholiq_p(coolantPressure_[i]), 0.49)) / 
                (pow(IF97::sigma97(IF97::Tsat97(coolantPressure_[i])), 0.5) * 
                pow(IF97::viscliq_p(coolantPressure_[i]), 0.29) * 
                pow(IF97::latentHeatVaporization_p(coolantPressure_[i]), 0.24) * 
                pow(IF97::rhovap_p(coolantPressure_[i]), 0.24)));

            hSB[i] = F * hc + S * hNB;
        }
        // Schrock-Grossman correlation (1962)
        else if (correlations_.lookup<word>("SaturatedBoiling") == "SchrockGrossman")
        {
            scalar tempQ(0.0);
            if (Q_[i] < 0.0)
            {
                tempQ = 0.0;
            }
            else
            {
                tempQ = Q_[i];
            }
            scalar recXtt(pow(tempQ / (1.0 - tempQ), 0.9) * 
                pow(IF97::rholiq_p(coolantPressure_[i]) / 
                IF97::rhovap_p(coolantPressure_[i]), 0.5) * 
                pow(IF97::viscvap_p(coolantPressure_[i]) / 
                IF97::viscliq_p(coolantPressure_[i]), 0.1));
            hSB[i] = (7400.0 * heatFlux_[i] / 
                IF97::latentHeatVaporization_p(coolantPressure_[i]) / 
                (massFlowRate_[i] / A_) + 1.11 * pow(recXtt, 0.66)) * hSP[i];
        }
        // Rohsenow correlation (1952)
        else if (correlations_.lookup<word>("SaturatedBoiling") == "Rohsenow")
        {
            hSB[i] = heatFlux_[i] / (cbrt(heatFlux_[i] / 
                (IF97::viscliq_p(coolantPressure_[i]) * 
                IF97::latentHeatVaporization_p(coolantPressure_[i]) * 
                pow(g_ * (IF97::rholiq_p(coolantPressure_[i]) - 
                IF97::rhovap_p(coolantPressure_[i])) / 
                IF97::sigma97(IF97::Tsat97(coolantPressure_[i])), 0.5))) *
                (0.013 * IF97::latentHeatVaporization_p(coolantPressure_[i]) * 
                pow(IF97::prandtlliq_p(coolantPressure_[i]), 1.0)) / 
                IF97::cpliq_p(coolantPressure_[i]) +
                IF97::Tsat97(coolantPressure_[i]) - T0_[i]);
        }
        else
        {
            FatalErrorInFunction()  
                << "Incorrectly specified saturated boiling HTC correlation name."  
                << abort(FatalError);
        }

        // Film boiling regime HTC
        // Groeneveld correlation (1973)
        if ((correlations_.lookup<word>("FilmBoiling") == "GroeneveldAnnular") || 
            (correlations_.lookup<word>("FilmBoiling") == "GroeneveldTubular"))
        {
            scalar tempQ(0.0);
            if (Q_[i] < 0.0)
            {
                tempQ = 0.0;
            }
            else
            {
                tempQ = Q_[i];
            }
            
            scalar GrA(0.0);
            scalar GrB(0.0);
            scalar GrC(0.0);
            scalar GrD(0.0);

            if (correlations_.lookup<word>("FilmBoiling") == "GroeneveldAnnular")
            {
                GrA = 0.052;
                GrB = 0.688;
                GrC = 1.26;
                GrD = -1.06;
            }
            else if (correlations_.lookup<word>("FilmBoiling") == "GroeneveldTubular")
            {
                GrA = 0.00109;
                GrB = 0.989;
                GrC = 1.41;
                GrD = -1.15;
            }
            else
            {
                FatalErrorInFunction()  
                    << "Incorrectly specified film boiling HTC correlation name (Groeneveld)."
                    << abort(FatalError);
            }

            scalar Y(max(1.0 - 0.1 * pow((1.0 - tempQ) * 
                (IF97::rholiq_p(coolantPressure_[i]) / 
                IF97::rhovap_p(coolantPressure_[i]) - 1.0), 0.4), 0.1));
            hFB[i] = GrA * IF97::tcondvap_p(coolantPressure_[i]) / DH_ * 
                pow(massFlowRate_[i] / A_ * DH_ / IF97::viscvap_p(coolantPressure_[i]) * 
                (tempQ + (1.0 - tempQ) * IF97::rhovap_p(coolantPressure_[i]) / 
                IF97::rholiq_p(coolantPressure_[i])), GrB) * 
                pow(Pr_[i], GrC) * pow(Y, GrD);
        }
        // Bishop-Sandberg-Tong correlation (1965)
        else if (correlations_.lookup<word>("FilmBoiling") == "BishopSandbergTong")
        {
            scalar rhoEq(X_[i] * IF97::rhovap_p(coolantPressure_[i]) +
                (1.0 - X_[i]) * IF97::rholiq_p(coolantPressure_[i]));
            hFB[i] = 0.0193 * IF97::tcond_Tp(Tfilm_[i], coolantPressure_[i]) / DH_ *
                pow(Re_[i], 0.8) * pow(Pr_[i], 1.23) * 
                pow(IF97::rhovap_p(coolantPressure_[i]) / rhoEq, 0.68) *
                pow(IF97::rhovap_p(coolantPressure_[i]) / 
                IF97::rholiq_p(coolantPressure_[i]), 0.068);
        }
        // Frederking correlation (1966)
        else if (correlations_.lookup<word>("FilmBoiling") == "Frederking")
        {
            scalar Tvap(0.5 * (Twall_[i] + IF97::Tsat97(coolantPressure_[i])));
            scalar hfg(IF97::latentHeatVaporization_p(coolantPressure_[i]) + 
                0.5 * IF97::cpvap_p(coolantPressure_[i]) * 
                (Twall_[i] - IF97::Tsat97(coolantPressure_[i])));
            scalar fsubc(1.0 + 0.1 * pow(IF97::rholiq_p(coolantPressure_[i]) / 
                IF97::rhovap_p(coolantPressure_[i]), 0.75) * 
                IF97::cpliq_p(coolantPressure_[i]) * 
                (IF97::Tsat97(coolantPressure_[i]) - T0_[i]) / 
                IF97::latentHeatVaporization_p(coolantPressure_[i]));
            hFB[i] = 0.2 * fsubc * cbrt(pow(IF97::tcond_Tp(Tvap, coolantPressure_[i]), 2.0) * 
                hfg * g_ * IF97::rhovap_p(coolantPressure_[i]) * 
                (IF97::rholiq_p(coolantPressure_[i]) - 
                IF97::rhovap_p(coolantPressure_[i])) / 
                IF97::viscvap_p(coolantPressure_[i]) /
                (Twall_[i] - IF97::Tsat97(coolantPressure_[i])));
        }
        // Sakurai correlation (1990)
        else if (correlations_.lookup<word>("FilmBoiling") == "Sakurai")
        {
            scalar Tvap(0.5 * (Twall_[i] + IF97::Tsat97(coolantPressure_[i])));
            scalar Tliq(0.5 * (T0_[i] + IF97::Tsat97(coolantPressure_[i])));
            scalar hcap(IF97::cpmass_Tp(Tvap, coolantPressure_[i]) * 
                (Twall_[i] - IF97::Tsat97(coolantPressure_[i])));
            scalar Lz(IF97::latentHeatVaporization_p(coolantPressure_[i]) + 0.5 * hcap);
            scalar Rz(sqrt((IF97::rhomass_Tp(Tvap, coolantPressure_[i]) * 
                IF97::visc_Tp(Tvap, coolantPressure_[i])) /
                (IF97::rhomass_Tp(Tliq, coolantPressure_[i]) * 
                IF97::visc_Tp(Tliq, coolantPressure_[i]))));
            scalar Prl(IF97::prandtl_Tp(Tliq, coolantPressure_[i]));
            scalar Prv(IF97::prandtl_Tp(Tvap, coolantPressure_[i]));
            scalar Sp(hcap / (Lz * Prv));
            scalar Sc(IF97::cpmass_Tp(Tliq, coolantPressure_[i]) * 
                (IF97::Tsat97(coolantPressure_[i]) - T0_[i]) / Lz);
            scalar Az(pow(Sc, 3.0) / 27.0 + Rz * Rz * Sp * Prl * Sc / 3.0 + 
                pow(Rz * Sp * Prl, 2.0) / 4.0);
            scalar Bz(-4.0 / 27.0 * Sc * Sc + 2.0 / 3.0 * Sp * Prl * Sc - 
                32.0 / 27.0 * Sp * Prl * Rz * Rz +
                0.25 * pow(Sp * Prl, 2.0) + 2.0 / 27.0 * pow(Sc, 3.0) / pow(Rz, 2.0));
            scalar Cz(0.5 * Rz * Rz * Sp * Prl);
            scalar Ez(cbrt(Az + Cz *sqrt(mag(Bz))) + cbrt(Az - Cz *sqrt(mag(Bz))) + Sc / 3.0);
            scalar lambdaTaylor(2.0 * constant::mathematical::pi * 
                sqrt(mag(IF97::sigma97(IF97::Tsat97(coolantPressure_[i])) /
                    (g_ * max(SMALL, IF97::rhomass_Tp(Tliq, coolantPressure_[i]) - 
                    IF97::rhomass_Tp(Tvap, coolantPressure_[i]))))));
            scalar Grv(g_ * (IF97::rhomass_Tp(Tliq, coolantPressure_[i]) - 
                IF97::rhomass_Tp(Tvap, coolantPressure_[i])) /
                (IF97::rhomass_Tp(Tvap, coolantPressure_[i]) * 
                pow(IF97::visc_Tp(Tvap, coolantPressure_[i]) / 
                IF97::rhomass_Tp(Tvap, coolantPressure_[i]), 2.0)) * 
                pow(lambdaTaylor, 3.0));
            scalar Mz((Grv * Prv * Lz * pow(Ez, 3.0)) / (pow(Rz * Prl * Sp, 2.0) * hcap *
                (1.0 + Ez / (Sp * Prl))));
            hFB[i] = 0.823 * IF97::tcond_Tp(Tvap, coolantPressure_[i]) / DH_ * pow(mag(Mz), 0.25);
        }
        else
        {
            FatalErrorInFunction()  
                << "Incorrectly specified film boiling HTC correlation name."  
                << abort(FatalError);
        }

        // Film boiling correction factor for RIA transients (equals to 1 by default)
        hFB[i] = hFB[i] * filmBoilingRIACorrection_;

        // Transition boiling regime HTC --> calculated after the evaluation of TCHF and qCHF_
    }

    // Evaluating ONB temperature
    scalarField TONB(this->size(), 600.0);

    // Basu correlation 
    // (N. Basu, G. R. Warrier, V. K. Dhir, Onset of Nucleate Boiling and Active Nucleation Site 
    // Density During Subcooled Flow Boiling, 2002)
    if (correlations_.lookup<word>("ONB") == "Basu")
    {
        forAll(zetas, i)
        {
            scalar phiAngle(0.663); // 0.995 for Zr-4, 0.663 for stainless steel
            scalar Fphi(1.0 - exp(-pow(phiAngle, 3.0) - 0.5 * phiAngle));
            scalar dTONB(2 * hSP[i] * IF97::Tsat97(coolantPressure_[i]) * 
                IF97::sigma97(IF97::Tsat97(coolantPressure_[i])) / 
                (pow(Fphi, 2.0) * IF97::rhovap_p(coolantPressure_[i]) * 
                IF97::latentHeatVaporization_p(coolantPressure_[i]) * coolantk_[i]));
            scalar dTsub(IF97::Tsat97(coolantPressure_[i]) - T0_[i]);

            TONB[i] = T0_[i] + 0.25 * pow((pow(dTONB, 0.5) + pow(dTONB + 4.0 * dTsub, 0.5)), 2.0);
        }
    }
    // Bergles-Rohsenow correlation 
    // (A. E. Bergles, W. M. Rohsenow, The Determination of Forced-Convection 
    // Surface-Boiling Heat Transfer, 1964)
    else if (correlations_.lookup<word>("ONB") == "BerglesRohsenow")
    {
        scalar oldValue(0.0);
        scalar fx(0.0);
        scalar dfx(0.0);
        forAll(zetas, i)
        {
            // Newton's method to solve for TONB
            do
            {
                fx = 15.6 * pow(coolantPressure_[i] / PsiToPa_, 1.156) * 
                    pow(mag(TONB[i] - IF97::Tsat97(coolantPressure_[i])) * 1.8, 2.3 / 
                    pow(coolantPressure_[i] / PsiToPa_, 0.0234)) * 
                    (BtuToJ_ / HrToSec_ / pow(FeetToMeters_, 2)) - 
                    hSP[i] * (TONB[i] - T0_[i]);

                dfx = 15.6 * pow(coolantPressure_[i] / PsiToPa_, 1.156) * 
                    pow(mag(TONB[i] - IF97::Tsat97(coolantPressure_[i])) * 1.8, 2.3 / 
                    pow(coolantPressure_[i] / PsiToPa_, 0.0234) - 1.0) * 
                    2.3 / pow(coolantPressure_[i] / PsiToPa_, 0.0234) * 
                    (BtuToJ_ / HrToSec_ / pow(FeetToMeters_, 2)) - hSP[i];

                oldValue = TONB[i];
                TONB[i] = oldValue - fx / dfx;
            }
            while(mag(TONB[i] - oldValue) / oldValue < 1e-5);
        }
    }
    // Bessiron correlation for NSRR RIA transients (stagnant flow, normal temperature and pressure) 
    // (V. Bessiron, Modelling of Clad-to-Coolant Heat Transfer for RIA Applications, 2007) 
    else if (correlations_.lookup<word>("ONB") == "Bessiron")
    {
        forAll(zetas, i)
        {
            TONB[i] = IF97::Tsat97(coolantPressure_[i]) + 20.0;
        }
    }
    else
    {
        FatalErrorInFunction()  
            << "Incorrectly specified ONB correlation name."  
            << abort(FatalError);
    }

    // Evaluate CHF heat flux
    scalarField TCHF(this->size(), 0.0);
    scalarField TTB(this->size(), 0.0);

    // EPRI-Columbia CHF correlation 
    // (D. G. Reddy, C. F. Fighetti, Parametric study of CHF data. Volume 2. A generalized 
    // subchannel CHF correlation for PWR and BWR fuel assemblies, NP-2609, Vol. 2, 1983)
    if (correlations_.lookup<word>("CHF") == "EPRI")
    {
        const scalar p1(0.5328);
        const scalar p2(0.1212);
        const scalar p3(1.6151);
        const scalar p4(1.4066);
        const scalar p5(-0.3040);
        const scalar p6(0.4843);
        const scalar p7(-0.3285);
        const scalar p8(-2.0749);
        scalar Qin((IF97::hmass_Tp(inletTemperature_, coolantPressure_[0]) - 
            IF97::hliq_p(coolantPressure_[0])) / 
            (IF97::hvap_p(coolantPressure_[0]) - IF97::hliq_p(coolantPressure_[0])));

        forAll(zetas, i)
        {
            scalar G(massFlowRate_[i] / A_ / 
                (1000000.0 * LbToKg_ / HrToSec_ / pow(FeetToMeters_, 2.0)));
            scalar reducedPressure(coolantPressure_[i] / 22064000.0);
            scalar Fa(pow(G, 0.1));
            scalar Fc(1.183 * Fa);
            scalar qzsum(0.0);
            scalar Y(0.0);
            for (int j = 0; j <= i; ++j)
            {
                qzsum = qzsum + heatFlux_[j] * faceHeights_[j];
            }
            if (heatFlux_[i] == 0.0)
            {
                Y = 0.0;
            }
            else
            {
                Y = qzsum / (heatFlux_[i] * zetas[i]);
            }
            scalar FAP(1 + (Y - 1.0) / (1.0 + G));
            scalar A(Fa * p1 * pow(reducedPressure, p2) * pow(G, p5 + p7 * reducedPressure));
            scalar C(Fc * FAP * p3 * pow(reducedPressure, p4) * pow(G, p6 + p8 * reducedPressure));
            if (heatFlux_[i] == 0.0)
            {
                qCHF_[i] = ROOTVGREAT;
            }
            else
            {
                qCHF_[i] = ((A - Qin) / (C + (Q_[i] - Qin) / (heatFlux_[i] / (1000000.0 * BtuToJ_ / 
                    HrToSec_ / pow(FeetToMeters_, 2.0))))) * (1000000.0 * BtuToJ_ / 
                    HrToSec_ / pow(FeetToMeters_, 2.0));
            }
        }
    }
    // Groeneveld CHF look-up table
    // (D. C. Groeneveld, The 2006 CHF look-up table, 
    // Nuclear Engineering and Design, Vol. 237, 2007)
    else if (correlations_.lookup<word>("CHF") == "lookuptableGroeneveld")
    {
        forAll(zetas, i)
        {
            // Interpolate between pressure tables
            FieldField<Field, scalar> pressureData
            (
                (*chfPressureTable_)(coolantPressure_[i] / 1000.0)
            );

            // Create 2D interpolate table (mass flow rate vs. quality)    
            scalarFieldInterpolateTable massFluxTable
            (
                chfMassFluxValues_, pressureData, chfMassFluxMethod_
            );

            // Interpolate between mass flux rows
            scalarField massFluxData(massFluxTable(massFlowRate_[i] / A_));

            // Create 1D quality interpolate table    
            scalarInterpolateTable qualityTable(chfQualityValues_, massFluxData, chfQualityMethod_);

            // Interpolate between quality columns, unit conversion to SI
            qCHF_[i] = (qualityTable(Q_[i]) * 1000.0);
            
            // Implement the correction (K1 for pipe diameter, K2 for rod bundle geometry)
            if (chfPitchToDiameter_ != 0.0)
            {
              scalar K2((2.0 * chfPitchToDiameter_ - 1.5) * exp(-cbrt(Q_[i]) / 2.0));
              qCHF_[i] = qCHF_[i] * K2;
            }
            else 
            {
              scalar K1(max(0.6, sqrt(0.008 / DH_)));
              qCHF_[i] = qCHF_[i] * K1;
            }
        }
    }
    // Zuber correlation 
    // (N. Zuber, Hydrodynamics Aspects of Boiling Heat Transfer, AECU-4439, 1959)
    else if (correlations_.lookup<word>("CHF") == "Zuber")
    {
        forAll(zetas, i)
        {
            qCHF_[i] = 0.131 * IF97::latentHeatVaporization_p(coolantPressure_[i]) * 
                pow(IF97::rhovap_p(coolantPressure_[i]), 0.5) * 
                pow(IF97::sigma97(IF97::Tsat97(coolantPressure_[i])) * g_ * 
                (IF97::rholiq_p(coolantPressure_[i]) - IF97::rhovap_p(coolantPressure_[i])), 0.25);
        }
    }
    // Modified Zuber correlation
    // (W. Liu, M. S. Kazimi, Modeling Cladding-Coolant Heat Transfer of
    // High-Burnup Fuel During RIA, 2006)
    else if (correlations_.lookup<word>("CHF") == "ModifiedZuber")
    {
        forAll(zetas, i)
        {
            scalar fsubc(1.0 + 0.1 * pow(IF97::rholiq_p(coolantPressure_[i]) / 
                IF97::rhovap_p(coolantPressure_[i]), 0.75) * IF97::cpliq_p(coolantPressure_[i]) * 
                (IF97::Tsat97(coolantPressure_[i]) - T0_[i]) / 
                IF97::latentHeatVaporization_p(coolantPressure_[i]));
            qCHF_[i] = 0.131 * fsubc * IF97::latentHeatVaporization_p(coolantPressure_[i]) * 
                pow(IF97::rhovap_p(coolantPressure_[i]), 0.5) *
                pow(IF97::sigma97(IF97::Tsat97(coolantPressure_[i])) * g_ * 
                (IF97::rholiq_p(coolantPressure_[i]) - IF97::rhovap_p(coolantPressure_[i])), 0.25);
        }
    }
    // GE new Hench-Levy CHF correlation 
    // (L. S. Tong, Boiling Crisis and Critical Heat Flux, TID-25887, 1972)
    else if (correlations_.lookup<word>("CHF") == "GE")
    {
        forAll(zetas, i)
        {
            scalar G(massFlowRate_[i] / A_ / (LbToKg_ / HrToSec_ / pow(FeetToMeters_, 2.0)));
            if (G < 500000.0)
            {
                qCHF_[i] = (840000.0 - 1000000.0 * Q_[i]) * 
                    (BtuToJ_ / HrToSec_ / pow(FeetToMeters_, 2.0));
            }
            else if (G < 750000.0)
            {
                qCHF_[i] = (800000.0 - 1000000.0 * Q_[i]) * 
                    (BtuToJ_ / HrToSec_ / pow(FeetToMeters_, 2.0));
            }
            else
            {
                qCHF_[i] = (800000.0 - 1000000.0 * Q_[i]) * 
                    (BtuToJ_ / HrToSec_ / pow(FeetToMeters_, 2.0));
                WarningInFunction() 
                    << "Mass flow rate is out of range for GE CHF correlation."  << endl;
            }

            if (heatFlux_[i] == 0.0)
            {
                qCHF_[i] = ROOTVGREAT;
            }
        }
    }
    // Bessiron correlation for NSRR RIA transients (stagnant flow, normal temperature and pressure)
    // (V. Bessiron, Modelling of Clad-to-Coolant Heat Transfer for RIA Applications, 2007) 
    else if (correlations_.lookup<word>("CHF") == "Bessiron")
    {
        forAll(zetas, i)
        {
            if (vaporizedThickness_[i] >= thicknessBessironRIA_)
            {
                qCHF_[i] = heatFlux_[i];
                vaporizedThickness_[i] = -1.0;
            }
        }
    }
    else
    {
        FatalErrorInFunction()  
            << "Incorrectly specified CHF correlation name."  
            << abort(FatalError);
    }

    // Evaluate CHF temperature
    forAll(zetas, i)
    {
        if (heatFlux_[i] == 0.0)
        {
            TCHF[i] = GREAT;
        }
        else if (correlations_.lookup<word>("CHF") == "Bessiron")
        {
            TCHF[i] = IF97::Tsat97(coolantPressure_[i]) + 20.0;
        }
        else
        {
            TCHF[i] = qCHF_[i] / hSB[i] + T0_[i];
        }

        // Transition boiling regime HTC
        // These correlations need to be calculated here as they require TCHF and qCHF
        // McDonough-Milich-King correlation (1961)
        if (correlations_.lookup<word>("TransitionBoiling") == "McDonoughMilichKing")
        {
            hTB[i] = heatFlux_[i] / ((qCHF_[i] - heatFlux_[i]) / 
                (4150.0 * exp(3970000.0 / coolantPressure_[i])) + TCHF[i] - 
                IF97::Tsat97(coolantPressure_[i]));
        }
        // Bessiron correlation for RIA transients (2007)
        else if (correlations_.lookup<word>("TransitionBoiling") == "Bessiron")
        {
            scalar SSH(hFB[i] * (Twall_[i] - IF97::Tsat97(coolantPressure_[i])));
            if (qCHF_[i] > 0.0)
            {
                TTB[i] = sqrt(mag(log(mag((heatFlux_[i] - SSH) / (qCHF_[i] - SSH))) / 
                    transitionBoilingRIACorrection_)) + TCHF[i];
            }
            else
            {
                TTB[i] = 0.0;
            }
            hTB[i] = max(heatFlux_[i] / (TTB[i] - IF97::Tsat97(coolantPressure_[i])), hFB[i]);
            if (Twall_[i] > IF97::Tsat97(coolantPressure_[i]) + 50.0)
            {
                hTB[i] = hFB[i];
            }
        }
        // Film boiling is assumed in the whole postCHF region
        else if (correlations_.lookup<word>("TransitionBoiling") == "None")
        {
            hTB[i] = hFB[i];
        }
        else
        {
            FatalErrorInFunction()  
                << "Incorrectly specified transition boiling HTC correlation name."  
                << abort(FatalError);
        }
    }


    // Evaluating Leidenfrost temperature and heat flux
    scalarField TLeidenfrost(this->size(), 0.0);
    scalarField qLeidenfrost(this->size(), 0.0);

    forAll(zetas, i)
    {
        // Corrected Groeneveld-Stewart correlation (1981)
        if (correlations_.lookup<word>("Leidenfrost") == "GroeneveldStewartCorrected")
        {
            TLeidenfrost[i] = 557.85 + 44.1 * coolantPressure_[i] / 1000000.0 -
                3.72 * pow(coolantPressure_[i] / 1000000.0, 2.0) -
                Q_[i] * 10000.0 / (2.82 + 1.22 * coolantPressure_[i] / 1000000.0);

            if (coolantPressure_[i] > 9e6)
            {
                TLeidenfrost[i] = (557.85 + 44.1 * 9000000.0 / 1000000.0 -
                    3.72 * pow(9000000.0 / 1000000.0, 2.0) -
                    Q_[i] * 10000.0 / (2.82 + 1.22 * 9000000.0 / 1000000.0) -
                    IF97::Tsat97(9000000.0)) * 
                    (22064000.0 - coolantPressure_[i]) / (22064000.0 - 9000000.0) +
                    IF97::Tsat97(coolantPressure_[i]);
            }

            qLeidenfrost[i] = hFB[i] * (TLeidenfrost[i] - IF97::Tsat97(coolantPressure_[i]));
        }
        // Groeneveld-Stewart correlation (1981)
        else if (correlations_.lookup<word>("Leidenfrost") == "GroeneveldStewart")
        {
            TLeidenfrost[i] = 557.85 + 44.1 * coolantPressure_[i] / 1000000.0 -
                3.72 * pow(coolantPressure_[i] / 1000000.0, 2.0);

            if (coolantPressure_[i] > 9e6)
            {
                TLeidenfrost[i] = (557.85 + 44.1 * 9000000.0 / 1000000.0 -
                    3.72 * pow(9000000.0 / 1000000.0, 2.0) - IF97::Tsat97(9000000.0)) * 
                    (22064000.0 - coolantPressure_[i]) / (22064000.0 - 9000000.0) +
                    IF97::Tsat97(coolantPressure_[i]);
            }

            qLeidenfrost[i] = hFB[i] * (TLeidenfrost[i] - IF97::Tsat97(coolantPressure_[i]));
        }
        // Ohnishi correlation (1984)
        else if (correlations_.lookup<word>("Leidenfrost") == "Ohnishi")
        {
            TLeidenfrost[i] = IF97::Tsat97(coolantPressure_[i]) + 350.0 + 
                5.1 * (IF97::Tsat97(coolantPressure_[i]) - T0_[i]);

            qLeidenfrost[i] = hFB[i] * (TLeidenfrost[i] - IF97::Tsat97(coolantPressure_[i]));
        }
        // Gotovskij correlation (1986)
        else if (correlations_.lookup<word>("Leidenfrost") == "Gotovskij")
        {
            TLeidenfrost[i] = IF97::Tsat97(coolantPressure_[i]) + 100.0 + 
                8.0 * (IF97::Tsat97(coolantPressure_[i]) - T0_[i]);
            qLeidenfrost[i] = hFB[i] * (TLeidenfrost[i] - IF97::Tsat97(coolantPressure_[i]));
        }
        // Zuber correlation (1959)
        else if (correlations_.lookup<word>("Leidenfrost") == "Zuber")
        {
            qLeidenfrost[i] = 0.176 * IF97::rhovap_p(coolantPressure_[i]) * 
                IF97::latentHeatVaporization_p(coolantPressure_[i]) *
                pow(g_ * IF97::sigma97(IF97::Tsat97(coolantPressure_[i])) * 
                    (IF97::rholiq_p(coolantPressure_[i]) - 
                    IF97::rhovap_p(coolantPressure_[i])) /
                    pow(IF97::rholiq_p(coolantPressure_[i]) + 
                    IF97::rhovap_p(coolantPressure_[i]), 2.0), 0.25);
            TLeidenfrost[i] = qLeidenfrost[i] / hFB[i] + IF97::Tsat97(coolantPressure_[i]);
        }
        // Bessiron correlation for RIA transients (2007)
        else if (correlations_.lookup<word>("Leidenfrost") == "Bessiron")
        {
            TLeidenfrost[i] = IF97::Tsat97(coolantPressure_[i]) + leidenfrostdTBessiron_;
            qLeidenfrost[i] = hFB[i] * (TLeidenfrost[i] - IF97::Tsat97(coolantPressure_[i]));
        }
        else
        {
            FatalErrorInFunction()  
                << "Incorrectly specified Leidenfrost point correlation name."  
                << abort(FatalError);
        }

        TLeidenfrost[i] = max(TLeidenfrost[i], IF97::Tsat97(coolantPressure_[i]) + 0.01);
    }

    forAll(zetas, i)
    {
        if (correlations_.lookup<word>("ONB") == "Bessiron" && correlationRegimeList_[i] == "saturatedBoiling" &&
         (vaporizedThickness_[i] >= 0.0 && vaporizedThickness_[i] < thicknessBessironRIA_))
        {
            Twall_[i] = TONB[i];
        }
    }
    
    // Evaluating heat transfer regime
    forAll(zetas, i)
    {
        // Single phase liquid
        if ((Q_[i] <= 0.0) && (Twall_[i] < TONB[i]))
        {
            h_[i] = hSP[i];
            correlationRegimeList_[i] = "singlePhaseLiquid";
        }
        // Single phase vapor
        else if (Q_[i] >= 1.0)
        {
            h_[i] = hSP[i];
            correlationRegimeList_[i] = "singlePhaseVapor";
        }
        // Subcooled boiling
        else if ((Q_[i] <= 0.0) && (Twall_[i] < TCHF[i]))
        {
            h_[i] = hSCB[i];
            correlationRegimeList_[i] = "subcooledBoiling";
        }
        // Saturated boiling
        else if (Twall_[i] < TCHF[i] || 
            (correlations_.lookup<word>("CHF") == "Bessiron" && vaporizedThickness_[i] >= 0.0 && 
            vaporizedThickness_[i] < thicknessBessironRIA_))
        {
            if (correlations_.lookup<word>("CHF") == "Bessiron" && 
                (this->patch().boundaryMesh().mesh().time().value() > timeCheck_))
            {
                if ((heatFlux_[i] - hSP[i] * (Twall_[i] - T0_[i])) > 0.0)
                {
                    vaporizedThickness_[i] = vaporizedThickness_[i] + 
                    (heatFlux_[i] - hSP[i] * (Twall_[i] - T0_[i])) / 
                    IF97::latentHeatVaporization_p(coolantPressure_[i]) /
                    IF97::rhomass_Tp(Tfilm_[i], coolantPressure_[i]) * 
                    this->patch().boundaryMesh().mesh().time().deltaT().value();
                }
            }

            h_[i] = hSB[i];
            correlationRegimeList_[i] = "saturatedBoiling";
        }
        // Post-CHF correlations
        else if ((Twall_[i] >= TCHF[i]) && (Q_[i] < 1.0))
        {
            // RIA post-CHF
            if (correlations_.lookup<word>("CHF") == "Bessiron" && heatFlux_[i] > qLeidenfrost[i])
            {
                h_[i] = hTB[i];
                correlationRegimeList_[i] = "RIApostCHF";
            }
            // Film boiling
            else if ((Twall_[i] >= TLeidenfrost[i] || qCHF_[i] <= heatFlux_[i]) 
                && heatFlux_[i] > qLeidenfrost[i])
            {
                h_[i] = hFB[i];
                correlationRegimeList_[i] = "filmBoiling";
            }
            // Transition boiling
            else if (Twall_[i] < TLeidenfrost[i] && heatFlux_[i] > qLeidenfrost[i])
            {
                h_[i] = hTB[i];
                correlationRegimeList_[i] = "transitionBoiling";
            }
            // Return to nucleate boiling
            else if (Twall_[i] < TLeidenfrost[i] || heatFlux_[i] <= qLeidenfrost[i])
            {
                h_[i] = hSB[i];
                correlationRegimeList_[i] = "saturatedBoiling";
            }
            else
            {
                FatalErrorInFunction()  
                    << "The simulation exited the known post-CHF fluid regimes."  
                    << abort(FatalError);
            }
        }
        else
        {
            FatalErrorInFunction()  
                << "The simulation exited the known fluid regimes."  
                << abort(FatalError);
        }
    }

    timeCheck_ = this->patch().boundaryMesh().mesh().time().value();

    // Relaxation of the HTC
    if (sum(oldH) > 0.0)
    {
        h_ = oldH * (1.0 - relaxHTC_) + h_ * relaxHTC_;
    }

    // Calculating alphas
    forAll(this->patch().Cf(),faceI)
    {
        alphaN[faceI] = h_[faceI];
    }

    alpha_ = alphaP / (alphaP + alphaN);

    if (correlations_.lookup<word>("ONB") == "Bessiron" && correlations_.lookup<word>("CHF") == "Bessiron")
    {
        scalarField alphaVap(this->size(), 0.0);
        alphaVap = (TONB - T0_) / (patchInternalField() - T0_ + SMALL);

        forAll(zetas, i)
        {
            if (correlationRegimeList_[i] == "saturatedBoiling" &&
            vaporizedThickness_[i] >= 0.0)
            {
                alpha_[i] = alphaVap[i];
                Twall_[i] = TONB[i];
            }
        }
    }

    Twall_ = alpha_*patchInternalField() + (1.0 - alpha_) * T0_;

    fvPatchField<scalar>::updateCoeffs();
}


void coolantChannelfvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<scalar>::operator=(
            alpha_*patchInternalField() + (1.0 - alpha_) * T0_
    );

    fvPatchField<scalar>::evaluate();
}


tmp<Field<scalar> >
coolantChannelfvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<scalarField>(new scalarField(alpha_));
}


tmp<Field<scalar> >
coolantChannelfvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return (1.0 - alpha_) * T0_;
}


tmp<Field<scalar> >
coolantChannelfvPatchScalarField::gradientInternalCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (alpha_ - 1.0) * deltaCoeffs;
}


tmp<Field<scalar> >
coolantChannelfvPatchScalarField::gradientBoundaryCoeffs() const
{
    const scalarField& deltaCoeffs = patch().deltaCoeffs();

    return (1.0 - alpha_) * deltaCoeffs * T0_;
}


void coolantChannelfvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
#ifdef OPENFOAMFOUNDATION

    if (!dict_.isNull())
    {
        writeEntry(os, "axialProfileDict", dict_);
    } 
    
    writeEntryIfDifferent<word>(os, "kappa", "k", kappaName_);
    writeEntry(os, "hydraulicDiameter", this->DH_);
    writeEntry(os, "enthalpyModel", this->enthalpyModel_);         
    writeEntry(os, "relaxHTC", this->relaxHTC_);
    writeEntry(os, "filmBoilingRIACorrection", this->filmBoilingRIACorrection_);
    writeEntry(os, "transitionBoilingRIACorrection", this->transitionBoilingRIACorrection_);
    writeEntry(os, "leidenfrostdTBessiron", this->leidenfrostdTBessiron_);
    writeEntry(os, "thicknessBessironRIA", this->thicknessBessironRIA_);
    writeEntry(os, "flowArea", this->A_);
    writeEntry(os, "flowDirection", this->flowDirection_);
    writeEntry(os, "correlations", correlations_);
    writeEntry(os, "correlationRegimeList", this->correlationRegimeList_);
    writeEntry(os, "T0", this->T0_);
    writeEntry(os, "heatFlux", this->heatFlux_);
    writeEntry(os, "value", this->Twall_);
    writeEntry(os, "coolantEnthalpy", this->coolantEnthalpy_);
    writeEntry(os, "voidFraction", this->X_);
    writeEntry(os, "steamQuality", this->Q_);
    writeEntry(os, "Re", this->Re_);
    writeEntry(os, "h", this->h_);
    if (correlations_.lookup<word>("CHF") == "Bessiron")
    {
        writeEntry(os, "vaporizedThickness", this->vaporizedThickness_);
        writeEntry(os, "qCHF", this->qCHF_);
    }
    if (coolantPressureList_.valid())
    {
        os.writeKeyword("coolantPressureList") << nl;
        os << token::BEGIN_BLOCK << nl;
        coolantPressureList_->write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        writeEntry<scalarField>(os, "coolantPressure", coolantPressure_);
    }
    if (inletTemperatureList_.valid())
    {
        os.writeKeyword("inletTemperatureList") << nl;
        os << token::BEGIN_BLOCK << nl;
        inletTemperatureList_->write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        writeEntry<scalar>(os, "inletTemperature", inletTemperature_);
    }
    if (massFlowRateList_.valid())
    {
        os.writeKeyword("massFlowRateList") << nl;
        os << token::BEGIN_BLOCK << nl;
        massFlowRateList_->write(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        writeEntry<scalarField>(os, "massFlowRate", massFlowRate_);
    }
#elif OPENFOAMESI

    if (!dict_.isNullDict())
    {
        os.writeEntry("axialProfileDict", dict_);      
    } 
    
    os.writeEntryIfDifferent<word>("kappa", "k", kappaName_);
    os.writeEntry("hydraulicDiameter", this->DH_);
    os.writeEntry("enthalpyModel", this->enthalpyModel_);
    os.writeEntry("relaxHTC", this->relaxHTC_);
    os.writeEntry("filmBoilingRIACorrection", this->filmBoilingRIACorrection_);
    os.writeEntry("transitionBoilingRIACorrection", this->transitionBoilingRIACorrection_);
    os.writeEntry("leidenfrostdTBessiron", this->leidenfrostdTBessiron_);
    os.writeEntry("thicknessBessironRIA", this->thicknessBessironRIA_);
    os.writeEntry("flowArea", this->A_);
    os.writeEntry("flowDirection", this->flowDirection_);
    os.writeEntry("correlations", correlations_);
    os.writeEntry("correlationRegimeList", this->correlationRegimeList_);
    os.writeEntry("T0", this->T0_);
    os.writeEntry("heatFlux"   , this->heatFlux_);
    os.writeEntry("value", this->Twall_);
    os.writeEntry("coolantEnthalpy", this->coolantEnthalpy_);   
    os.writeEntry("voidFraction", this->X_);
    os.writeEntry("steamQuality", this->Q_);
    os.writeEntry("Re", this->Re_);
    os.writeEntry("h", this->h_);
    if (correlations_.lookup<word>("CHF") == "Bessiron")
    {
        os.writeEntry("vaporizedThickness", this->vaporizedThickness_);
        os.writeEntry("qCHF", this->qCHF_);
    }
    if (coolantPressureList_.valid())
    {
        os.writeKeyword("coolantPressureList") << nl;
        os << token::BEGIN_BLOCK << nl;
        coolantPressureList_->writeData(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        coolantPressure_.writeEntry("coolantPressure", os);
    }
    if (inletTemperatureList_.valid())
    {
        os.writeKeyword("inletTemperatureList") << nl;
        os << token::BEGIN_BLOCK << nl;
        inletTemperatureList_->writeData(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        os.writeEntry("inletTemperature", this->inletTemperature_);
    }
    if (massFlowRateList_.valid())
    {
        os.writeKeyword("massFlowRateList") << nl;
        os << token::BEGIN_BLOCK << nl;
        massFlowRateList_->writeData(os);
        os << token::END_BLOCK << nl;
    }
    else
    {
        massFlowRate_.writeEntry("massFlowRate", os);
    }
#endif  
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, coolantChannelfvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
