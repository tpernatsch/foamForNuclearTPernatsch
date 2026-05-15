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

#include "coolantChannelModel.H"
#include "IF97.H"
#include "Pstream.H"
#include "error.H"
#include "mathematicalConstants.H"
#include "ListOps.H"

namespace Foam
{

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void coolantChannelModel::initializeAxialProfile(const dictionary& dict)
{
    const dictionary& axialDict = dict.subDict("axialProfileDict");

    T0ProfileData_ = FieldField<Field, scalar>
    (
        PtrList<scalarField>(axialDict.lookup("T0Data"), scalarFieldFieldINew())
    );

    zMethod_ = interpolateTableBase::interpolationMethod
    (
        interpolateTableBase::interpolationMethodNames_
        [
            axialDict.lookupOrDefault<word>("axialInterpolationMethod", "linear")
        ]
    );

    tMethod_ = interpolateTableBase::interpolationMethod
    (
        interpolateTableBase::interpolationMethodNames_
        [
            axialDict.lookupOrDefault<word>("timeInterpolationMethod", "linear")
        ]
    );

    zValues_ = scalarField(axialDict.lookup("axialLocations"));
    timeValues_ = scalarField(axialDict.lookup("timePoints"));
    axialProfileDict_ = dictionary(axialDict);

    initializeAxialProfileTable();
}


void coolantChannelModel::initializeAxialProfileTable()
{
    if (enthalpyModel_ || axialProfileDict_.isNullDict())
    {
        return;
    }

#ifdef OPENFOAMFOUNDATION
    T0Table_.set
    (
        new scalarFieldInterpolateTable
        (
            timeValues_,
            T0ProfileData_,
            tMethod_
        )
    );
#elif OPENFOAMESI
    T0Table_.reset
    (
        new scalarFieldInterpolateTable
        (
            timeValues_,
            T0ProfileData_,
            tMethod_
        )
    );
#endif
}


void coolantChannelModel::validate() const
{
    if (A_ <= SMALL)
    {
        FatalErrorInFunction
            << "flowArea must be > 0. Current value = " << A_
            << abort(FatalError);
    }

    if (DH_ <= SMALL)
    {
        FatalErrorInFunction
            << "hydraulicDiameter must be > 0. Current value = " << DH_
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

coolantChannelModel::coolantChannelModel(const label nFaces)
:
    A_(-1.0),
    flowDirection_(1.0),
    DH_(-1.0),
    coolantPressure_(nFaces, 0.0),
    coolantPressureList_(),
    massFlowRate_(nFaces, -1.0),
    massFlowRateList_(),
    inletTemperature_(-1.0),
    inletTemperatureList_(),
    enthalpyModel_(true),
    inletTemperatureRequired_(false),
    positiveMassFlowRequired_(true),
    zValues_(),
    timeValues_(),
    T0ProfileData_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    T0Table_(),
    axialProfileDict_(dictionary::null),
    T0_(nFaces, 0.0),
    Q_(nFaces, 0.0),
    X_(nFaces, 0.0),
    Re_(nFaces, 0.0),
    Pr_(nFaces, 0.0),
    Tfilm_(nFaces, 0.0),
    coolantk_(nFaces, 0.0),
    coolantCp_(nFaces, 0.0),
    coolantEnthalpy_(nFaces, 0.0),
    rho_(nFaces, 0.0),
    rhoOld_(nFaces, 0.0),
    hOld_(nFaces, 0.0),
    velocity_(nFaces, 0.0),
    maxPicardIter_(10),
    picardTol_(1.0),
    timeCheckHEM_(-GREAT)
{}


coolantChannelModel::coolantChannelModel
(
    const dictionary& dict,
    const label nFaces,
    const bool inletTemperatureRequired,
    const bool enthalpyModelDefault,
    const bool positiveMassFlowRequired
)
:
    A_(readScalar(dict.lookup("flowArea"))),
    flowDirection_(dict.lookupOrDefault<scalar>("flowDirection", 1.0)),
    DH_(readScalar(dict.lookup("hydraulicDiameter"))),
    coolantPressure_(nFaces, 0.0),
    coolantPressureList_(),
    massFlowRate_(nFaces, -1.0),
    massFlowRateList_(),
    inletTemperature_(-1.0),
    inletTemperatureList_(),
#ifdef OPENFOAMFOUNDATION
    enthalpyModel_
    (
        dict.lookupOrDefault<bool>("enthalpyModel", enthalpyModelDefault)
    ),
#elif OPENFOAMESI
    enthalpyModel_
    (
        dict.getOrDefault<bool>("enthalpyModel", enthalpyModelDefault)
    ),
#endif
    inletTemperatureRequired_(inletTemperatureRequired),
    positiveMassFlowRequired_(positiveMassFlowRequired),
    zValues_(),
    timeValues_(),
    T0ProfileData_(),
    zMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    tMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    T0Table_(),
    axialProfileDict_(dictionary::null),
    T0_(nFaces, 0.0),
    Q_(nFaces, 0.0),
    X_(nFaces, 0.0),
    Re_(nFaces, 0.0),
    Pr_(nFaces, 0.0),
    Tfilm_(nFaces, 0.0),
    coolantk_(nFaces, 0.0),
    coolantCp_(nFaces, 0.0),
    coolantEnthalpy_(nFaces, 0.0),
    rho_(nFaces, 0.0),
    rhoOld_(nFaces, 0.0),
    hOld_(nFaces, 0.0),
    velocity_(nFaces, 0.0),
    maxPicardIter_(dict.lookupOrDefault<label>("maxPicardIter", 10)),
    picardTol_(dict.lookupOrDefault<scalar>("picardTol", 1.0)),
    timeCheckHEM_(-GREAT)
{
    if (dict.found("coolantPressureList"))
    {
        coolantPressureList_.reset
        (
            new scalarTable
            (
                "coolantPressureList",
                dict.subDict("coolantPressureList")
            )
        );
    }
    else if (dict.found("coolantPressure"))
    {
        coolantPressure_ = scalarField("coolantPressure", dict, nFaces);
    }
    else
    {
        FatalErrorInFunction
            << "Missing coolant input. Provide either \"coolantPressureList\" "
            << "or \"coolantPressure\"."
            << abort(FatalError);
    }

    if (dict.found("massFlowRateList"))
    {
        massFlowRateList_.reset
        (
            new scalarTable
            (
                "massFlowRateList",
                dict.subDict("massFlowRateList")
            )
        );
    }
    else if (dict.found("massFlowRate"))
    {
        massFlowRate_ = scalarField("massFlowRate", dict, nFaces);
    }
    else
    {
        FatalErrorInFunction
            << "Missing coolant input. Provide either \"massFlowRateList\" "
            << "or \"massFlowRate\"."
            << abort(FatalError);
    }

    if (enthalpyModel_ || inletTemperatureRequired_)
    {
        if (dict.found("inletTemperatureList"))
        {
            inletTemperatureList_.reset
            (
                new scalarTable
                (
                    "inletTemperatureList",
                    dict.subDict("inletTemperatureList")
                )
            );
        }
        else if (dict.found("inletTemperature"))
        {
            inletTemperature_ = readScalar(dict.lookup("inletTemperature"));
        }
        else
        {
            FatalErrorInFunction
                << "enthalpyModel=true or chf=EPRI requires either "
                << "\"inletTemperatureList\" or \"inletTemperature\"."
                << abort(FatalError);
        }
    }

    if (!enthalpyModel_)
    {
        if (!dict.found("axialProfileDict"))
        {
            FatalErrorInFunction
                << "enthalpyModel=false requires \"axialProfileDict\"."
                << abort(FatalError);
        }

        initializeAxialProfile(dict);
    }

    // Seed coolant state from stored values (restart) or inlet conditions so
    // that the first write() before any updateCoeffs() emits meaningful values.
    if (dict.found("T0"))
    {
        T0_ = scalarField("T0", dict, nFaces);
    }
    else if (inletTemperature_ > 0)
    {
        T0_ = inletTemperature_;
    }

    if (enthalpyModel_ && inletTemperature_ > 0)
    {
        if (dict.found("coolantEnthalpy"))
        {
            coolantEnthalpy_ = scalarField("coolantEnthalpy", dict, nFaces);
            hOld_ = coolantEnthalpy_;
        }
        else
        {
            // Use uniform inlet enthalpy as first-write placeholder.
            // First-call detection in updateBulkState will re-initialise
            // properly on the first solve.
            const scalar hInlet =
                IF97::hmass_Tp(inletTemperature_, coolantPressure_[0]);
            coolantEnthalpy_ = hInlet;
            hOld_            = hInlet;
        }

        if (dict.found("rho"))
        {
            rho_    = scalarField("rho", dict, nFaces);
            rhoOld_ = rho_;
        }
        else
        {
            const scalar rhoInlet =
                IF97::rhomass_phmass(coolantPressure_[0], coolantEnthalpy_[0]);
            rho_    = rhoInlet;
            rhoOld_ = rhoInlet;
        }
    }

    validate();
}


coolantChannelModel::coolantChannelModel(const coolantChannelModel& model)
:
    A_(model.A_),
    flowDirection_(model.flowDirection_),
    DH_(model.DH_),
    coolantPressure_(model.coolantPressure_),
    coolantPressureList_(),
    massFlowRate_(model.massFlowRate_),
    massFlowRateList_(),
    inletTemperature_(model.inletTemperature_),
    inletTemperatureList_(),
    enthalpyModel_(model.enthalpyModel_),
    inletTemperatureRequired_(model.inletTemperatureRequired_),
    positiveMassFlowRequired_(model.positiveMassFlowRequired_),
    zValues_(model.zValues_),
    timeValues_(model.timeValues_),
    T0ProfileData_(model.T0ProfileData_),
    zMethod_(model.zMethod_),
    tMethod_(model.tMethod_),
    T0Table_(),
    axialProfileDict_(model.axialProfileDict_),
    T0_(model.T0_),
    Q_(model.Q_),
    X_(model.X_),
    Re_(model.Re_),
    Pr_(model.Pr_),
    Tfilm_(model.Tfilm_),
    coolantk_(model.coolantk_),
    coolantCp_(model.coolantCp_),
    coolantEnthalpy_(model.coolantEnthalpy_),
    rho_(model.rho_),
    rhoOld_(model.rhoOld_),
    hOld_(model.hOld_),
    velocity_(model.velocity_),
    maxPicardIter_(model.maxPicardIter_),
    picardTol_(model.picardTol_),
    timeCheckHEM_(model.timeCheckHEM_)
{
    initializeAxialProfileTable();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coolantChannelModel::updateBoundaryInputs(const scalar userTime)
{
    if (coolantPressureList_.valid())
    {
        coolantPressure_ = coolantPressureList_->value(userTime);
    }

    if (massFlowRateList_.valid())
    {
        massFlowRate_ = massFlowRateList_->value(userTime);
    }

    if (inletTemperatureList_.valid())
    {
        inletTemperature_ = inletTemperatureList_->value(userTime);
    }

    forAll(massFlowRate_, faceI)
    {
        if (positiveMassFlowRequired_ && massFlowRate_[faceI] <= SMALL)
        {
            FatalErrorInFunction
                << "massFlowRate must be > 0 for all faces. Found "
                << massFlowRate_[faceI] << " on face " << faceI
                << abort(FatalError);
        }
    }

    forAll(coolantPressure_, faceI)
    {
        if (coolantPressure_[faceI] <= SMALL)
        {
            FatalErrorInFunction
                << "coolantPressure must be > 0 for all faces. Found "
                << coolantPressure_[faceI] << " on face " << faceI
                << abort(FatalError);
        }
    }
}


void coolantChannelModel::updateBulkState
(
    const scalar userTime,
    const scalarField& zetas,
    const scalarField& radii,
    const scalarField& faceHeights,
    const scalarField& heatFlux,
    const scalar deltaT
)
{
    if (enthalpyModel_)
    {
        const scalar g = 9.80665;

        // --- Freeze old-time state at the start of each new time step ---
        if (userTime > timeCheckHEM_)
        {
            hOld_ = coolantEnthalpy_;
            rhoOld_ = rho_;
            timeCheckHEM_ = userTime;
        }

        // --- Gather all face data across processors ---
        // The global array is ordered proc0, proc1, ..., so local face faceI
        // maps to global index (procOffset + faceI).

        scalarField allZetas, allFaceHeights, allRadii;
        scalarField allHeatFlux, allPressure, allHOld, allRhoOld, allH;
        label procOffset = 0;

        if (Pstream::parRun())
        {
            List<scalarField> gZ(Pstream::nProcs(), zetas);
            List<scalarField> gFH(Pstream::nProcs(), faceHeights);
            List<scalarField> gR(Pstream::nProcs(), radii);
            List<scalarField> gQ(Pstream::nProcs(), heatFlux);
            List<scalarField> gP(Pstream::nProcs(), coolantPressure_);
            List<scalarField> gHOld(Pstream::nProcs(), hOld_);
            List<scalarField> gRhoOld(Pstream::nProcs(), rhoOld_);
            List<scalarField> gH(Pstream::nProcs(), coolantEnthalpy_);

            Pstream::gatherList(gZ);      Pstream::scatterList(gZ);
            Pstream::gatherList(gFH);     Pstream::scatterList(gFH);
            Pstream::gatherList(gR);      Pstream::scatterList(gR);
            Pstream::gatherList(gQ);      Pstream::scatterList(gQ);
            Pstream::gatherList(gP);      Pstream::scatterList(gP);
            Pstream::gatherList(gHOld);   Pstream::scatterList(gHOld);
            Pstream::gatherList(gRhoOld); Pstream::scatterList(gRhoOld);
            Pstream::gatherList(gH);      Pstream::scatterList(gH);

            for (label procI = 0; procI < Pstream::myProcNo(); ++procI)
            {
                procOffset += gZ[procI].size();
            }

            forAll(gZ, procI)
            {
                allZetas.append(gZ[procI]);
                allFaceHeights.append(gFH[procI]);
                allRadii.append(gR[procI]);
                allHeatFlux.append(gQ[procI]);
                allPressure.append(gP[procI]);
                allHOld.append(gHOld[procI]);
                allRhoOld.append(gRhoOld[procI]);
                allH.append(gH[procI]);
            }
        }
        else
        {
            allZetas       = zetas;
            allFaceHeights = faceHeights;
            allRadii       = radii;
            allHeatFlux    = heatFlux;
            allPressure    = coolantPressure_;
            allHOld        = hOld_;
            allRhoOld      = rhoOld_;
            allH           = coolantEnthalpy_;
            procOffset     = 0;
        }

        const label nAll = allZetas.size();

        // --- On first call: initialise from inlet conditions ---
        if (max(allRhoOld) < SMALL)
        {
            forAll(allRhoOld, i)
            {
                const scalar hInlet =
                    IF97::hmass_Tp(inletTemperature_, allPressure[i]);
                allHOld[i]   = hInlet;
                allH[i]      = hInlet;
                allRhoOld[i] = IF97::rhomass_phmass(allPressure[i], hInlet);
            }
        }

        // --- Sort faces by axial position in flow direction ---
        scalarField sortZ(nAll);
        forAll(sortZ, i) { sortZ[i] = allZetas[i]*flowDirection_; }
        const labelList sortOrder = sortedOrder(sortZ);

        // --- Inlet mass flow rate ---
        // massFlowRate_ was loaded uniformly by updateBoundaryInputs;
        // take element 0 as the inlet value.
        const scalar mdotInlet = massFlowRate_[0];

        // Results indexed by global face index
        scalarField marchH(nAll, 0.0);
        scalarField marchRho(nAll, 0.0);
        scalarField marchV(nAll, 0.0);
        scalarField marchMdot(nAll, 0.0);

        // --- Sequential Picard march from inlet to outlet ---
        // Inlet total enthalpy: h_inlet + gz_inlet (no KE at inlet assumed)
        const label iFirst = sortOrder[0];
        scalar mdotIn = mdotInlet;
        scalar H_up   = IF97::hmass_Tp(inletTemperature_, allPressure[iFirst])
                      + g*allZetas[iFirst];

        forAll(sortOrder, k)
        {
            const label i = sortOrder[k];

            const scalar p_i = allPressure[i];
            const scalar Dz  = allFaceHeights[i];
            const scalar P_i = 2.0*constant::mathematical::pi*allRadii[i];
            const scalar z_i = allZetas[i];
            const scalar q_i = allHeatFlux[i];
            const scalar rhoN = allRhoOld[i];
            const scalar hN   = allHOld[i];

            // Backward-Euler coefficient: ρ^n * A * Δz / Δt
            const scalar coeff = rhoN*A_*Dz / deltaT;

            // Wall heat source
            const scalar src = q_i*P_i*Dz;

            // Picard iteration (warm start from previous march result)
            scalar h_k = allH[i];

            for (label iter = 0; iter < maxPicardIter_; ++iter)
            {
                const scalar rho_k = IF97::rhomass_phmass(p_i, h_k);

                // Continuity: outlet mdot from density change
                const scalar mdotOut =
                    mdotIn - (rho_k - rhoN)*A_*Dz/deltaT;

                const scalar v_k = mdotOut / max(rho_k*A_, SMALL);

                // Outlet KE + PE (upwind: outlet face value = cell centre)
                const scalar KEpe = 0.5*sqr(v_k) + g*z_i;

                // Solve energy equation linearly for h^{n+1}
                const scalar denom = max(coeff + mdotOut, SMALL);
                const scalar h_new =
                    (coeff*hN + mdotIn*H_up + src - mdotOut*KEpe)
                  / denom;

                const scalar residual = mag(h_new - h_k);
                h_k = h_new;

                if (residual < picardTol_)
                {
                    break;
                }
            }

            // Commit converged state for face i
            const scalar rho_f   = IF97::rhomass_phmass(p_i, h_k);
            const scalar mdotOut = mdotIn - (rho_f - rhoN)*A_*Dz/deltaT;
            const scalar v_f     = mdotOut / max(rho_f*A_, SMALL);

            marchH[i]    = h_k;
            marchRho[i]  = rho_f;
            marchV[i]    = v_f;
            marchMdot[i] = mdotOut;

            // Pass downstream
            mdotIn = mdotOut;
            H_up   = h_k + 0.5*sqr(v_f) + g*z_i;
        }

        // --- Write results to local fields ---
        // Local face faceI corresponds to global index (procOffset + faceI)
        forAll(zetas, faceI)
        {
            const label gi = procOffset + faceI;
            coolantEnthalpy_[faceI] = marchH[gi];
            rho_[faceI]             = marchRho[gi];
            velocity_[faceI]        = marchV[gi];
            massFlowRate_[faceI]    = marchMdot[gi];
            T0_[faceI] = IF97::T_phmass(allPressure[gi], marchH[gi]);
        }
    }
    else
    {
        scalarField T0zData(T0Table_()(userTime));

        scalarInterpolateTable T0zTable
        (
            zValues_,
            T0zData,
            zMethod_,
            interpolateTableBase::outOfBoundsMethod::FIXED
        );

        forAll(zetas, faceI)
        {
            T0_[faceI] = T0zTable(zetas[faceI]);
            coolantEnthalpy_[faceI] =
                IF97::hmass_Tp(T0_[faceI], coolantPressure_[faceI]);
        }
    }
}


void coolantChannelModel::updateDerivedState(const scalarField& Twall)
{
    forAll(T0_, i)
    {
        Q_[i] =
            (coolantEnthalpy_[i] - IF97::hliq_p(coolantPressure_[i]))
          /
            (
                IF97::hvap_p(coolantPressure_[i])
              - IF97::hliq_p(coolantPressure_[i])
            );

        Tfilm_[i] = 0.5*(T0_[i] + Twall[i]);

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
            X_[i] =
                1.0
              /
                (
                    1.0
                  + (1.0 - Q_[i])/Q_[i]
                  * IF97::rhovap_p(coolantPressure_[i])
                  / IF97::rholiq_p(coolantPressure_[i])
                );
        }

        Re_[i] =
            massFlowRate_[i]*DH_
          / IF97::visc_Tp(Tfilm_[i], coolantPressure_[i])
          / A_;

        Pr_[i] = IF97::prandtl_Tp(Tfilm_[i], coolantPressure_[i]);
        coolantk_[i] = IF97::tcond_Tp(Tfilm_[i], coolantPressure_[i]);
        coolantCp_[i] = IF97::cpmass_Tp(Tfilm_[i], coolantPressure_[i]);
    }
}


void coolantChannelModel::write(Ostream& os) const
{
    if (!enthalpyModel_ && !axialProfileDict_.isNullDict())
    {
        os.writeEntry("axialProfileDict", axialProfileDict_);
    }

    os.writeEntry("hydraulicDiameter", DH_);
    os.writeEntry("flowArea", A_);
    os.writeEntry("flowDirection", flowDirection_);
    os.writeEntry("enthalpyModel", enthalpyModel_);

    T0_.writeEntry("T0", os);
    coolantEnthalpy_.writeEntry("coolantEnthalpy", os);
    X_.writeEntry("voidFraction", os);
    Q_.writeEntry("steamQuality", os);
    Re_.writeEntry("Re", os);
    Pr_.writeEntry("Pr", os);

    if (enthalpyModel_)
    {
        os.writeEntry("maxPicardIter", maxPicardIter_);
        os.writeEntry("picardTol", picardTol_);
        rho_.writeEntry("rho", os);
        velocity_.writeEntry("velocity", os);
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

    if (enthalpyModel_ || inletTemperatureRequired_)
    {
        if (inletTemperatureList_.valid())
        {
            os.writeKeyword("inletTemperatureList") << nl;
            os << token::BEGIN_BLOCK << nl;
            inletTemperatureList_->writeData(os);
            os << token::END_BLOCK << nl;
        }
        else
        {
            os.writeEntry("inletTemperature", inletTemperature_);
        }
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
}

} // End namespace Foam

// ************************************************************************* //
