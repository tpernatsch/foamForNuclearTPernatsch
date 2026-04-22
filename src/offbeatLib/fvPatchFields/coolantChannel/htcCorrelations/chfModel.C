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

#include "chfModel.H"
#include "IF97.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "IOobject.H"
#include "scalarFieldFieldINew.H"
#include "error.H"

#ifdef OPENFOAMFOUNDATION
    // nothing extra
#elif OPENFOAMESI
    // nothing extra
#endif

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

const scalar chfModel::PsiToPa_     = 6894.757;
const scalar chfModel::BtuToJ_      = 1055.06;
const scalar chfModel::FeetToMeters_ = 0.3048;
const scalar chfModel::LbToKg_      = 0.45359237;
const scalar chfModel::HrToSec_     = 3600.0;
const scalar chfModel::g_           = 9.80665;


chfModel::Model
chfModel::parse(const word& name)
{
    if (name == "EPRI")           return Model::EPRI;
    if (name == "GroeneveldTable") return Model::GroeneveldTable;
    if (name == "Zuber")          return Model::Zuber;
    if (name == "ModifiedZuber")  return Model::ModifiedZuber;
    if (name == "GE")             return Model::GE;

    FatalErrorInFunction
        << "Unknown CHF correlation \"" << name << "\"." << nl
        << "Valid options: EPRI, GroeneveldTable, Zuber, ModifiedZuber, GE"
        << abort(FatalError);

    return Model::EPRI;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

chfModel::chfModel(const word& name)
:
    model_(parse(name)),
    Qin_(0.0),
    FAP_(),
    chfData_(),
    chfPitchToDiameter_(0.0),
    chfPressureValues_(),
    chfMassFluxValues_(),
    chfQualityValues_(),
    chfPressureMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    chfMassFluxMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    chfQualityMethod_(interpolateTableBase::interpolationMethodNames_["linear"]),
    chfPressureTable_()
{}


chfModel::chfModel(const chfModel& other)
:
    model_(other.model_),
    Qin_(other.Qin_),
    FAP_(other.FAP_),
    chfData_(other.chfData_),
    chfPitchToDiameter_(other.chfPitchToDiameter_),
    chfPressureValues_(other.chfPressureValues_),
    chfMassFluxValues_(other.chfMassFluxValues_),
    chfQualityValues_(other.chfQualityValues_),
    chfPressureMethod_(other.chfPressureMethod_),
    chfMassFluxMethod_(other.chfMassFluxMethod_),
    chfQualityMethod_(other.chfQualityMethod_),
    chfPressureTable_()
{
    // Rebuild the interpolation wrapper from the deep-copied data.
    // We never copy chfPressureTable_ directly because autoPtr does not
    // support copying (the original would be left empty).
    if (other.chfPressureTable_.valid())
    {
        #ifdef OPENFOAMFOUNDATION
        chfPressureTable_.set
        (
            new scalarFieldFieldInterpolateTable
            (
                chfPressureValues_,
                chfData_,
                chfPressureMethod_
            )
        );
        #elif OPENFOAMESI
        chfPressureTable_.reset
        (
            new scalarFieldFieldInterpolateTable
            (
                chfPressureValues_,
                chfData_,
                chfPressureMethod_
            )
        );
        #endif
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

scalar chfModel::computeInletQuality
(
    const scalar inletTemperature,
    const scalar pInlet
)
{
    return
        (
            IF97::hmass_Tp(inletTemperature, pInlet)
          - IF97::hliq_p(pInlet)
        )
      / (IF97::hvap_p(pInlet) - IF97::hliq_p(pInlet));
}


scalarField chfModel::computeShapeFactors
(
    const scalarField& zetas,
    const scalarField& faceHeights,
    const scalarField& heatFlux,
    const scalarField& massFlowRate,
    const scalar flowArea,
    const scalar flowDirection
)
{
    const label nFaces = zetas.size();

    scalarField FAP(nFaces, 1.0);

    // G in Mlb/(h·ft²) — uniform across faces (use first face value)
    const scalar GMlb =
        massFlowRate[0]/flowArea
      / (1e6 * LbToKg_/HrToSec_/pow(FeetToMeters_, 2.0));

    forAll(zetas, faceI)
    {
        const scalar zI = zetas[faceI]*sign(flowDirection);

        // Cumulative heat flux integral from inlet to face faceI
        scalar qzsum = 0.0;
        forAll(zetas, j)
        {
            if (zetas[j]*sign(flowDirection) <= zI)
            {
                qzsum += heatFlux[j]*faceHeights[j];
            }
        }

        const scalar zAbs = mag(zetas[faceI]);
        scalar Y = 0.0;

        if (mag(heatFlux[faceI]) > SMALL && zAbs > SMALL)
        {
            Y = qzsum/(heatFlux[faceI]*zAbs);
        }

        FAP[faceI] = 1.0 + (Y - 1.0)/(1.0 + GMlb);
    }

    return FAP;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void chfModel::precompute
(
    const scalarField& zetas,
    const scalarField& faceHeights,
    const scalarField& heatFlux,
    const scalarField& massFlowRate,
    const scalarField& coolantPressure,
    const scalar inletTemperature,
    const scalar flowArea,
    const scalar flowDirection
)
{
    if (model_ != Model::EPRI)
    {
        return;
    }

    // Find the inlet face (minimum axial position in the flow direction)
    label inletFace = 0;
    scalar zInlet = zetas[0]*sign(flowDirection);

    for (label k = 1; k < zetas.size(); ++k)
    {
        const scalar zk = zetas[k]*sign(flowDirection);
        if (zk < zInlet)
        {
            zInlet = zk;
            inletFace = k;
        }
    }

    Qin_ = computeInletQuality(inletTemperature, coolantPressure[inletFace]);

    FAP_ = computeShapeFactors
    (
        zetas,
        faceHeights,
        heatFlux,
        massFlowRate,
        flowArea,
        flowDirection
    );
}


void chfModel::initializeTable
(
    const dictionary& bcDict,
    const fvMesh& mesh
)
{
    IOdictionary chfDict
    (
        IOobject
        (
            "chfTableDict",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    chfData_ = PtrList<FieldField<Field, scalar>>
    (
        PtrList<FieldField<Field, scalar>>
        (
            chfDict.lookup("data"),
            PtrListScalarFieldFieldINew()
        )
    );

    chfPitchToDiameter_ =
        chfDict.lookupOrDefault<scalar>("PitchToDiameter", 0.0);

    chfPressureValues_ = scalarField(chfDict.lookup("pressureValues"));
    chfMassFluxValues_ = scalarField(chfDict.lookup("massFluxValues"));
    chfQualityValues_  = scalarField(chfDict.lookup("qualityValues"));

    chfPressureMethod_ = interpolateTableBase::interpolationMethod
    (
        interpolateTableBase::interpolationMethodNames_
        [
            bcDict.lookupOrDefault<word>("pressureInterpolationMethod", "linear")
        ]
    );

    chfMassFluxMethod_ = interpolateTableBase::interpolationMethod
    (
        interpolateTableBase::interpolationMethodNames_
        [
            bcDict.lookupOrDefault<word>("massFluxInterpolationMethod", "linear")
        ]
    );

    chfQualityMethod_ = interpolateTableBase::interpolationMethod
    (
        interpolateTableBase::interpolationMethodNames_
        [
            bcDict.lookupOrDefault<word>("qualityInterpolationMethod", "linear")
        ]
    );

    #ifdef OPENFOAMFOUNDATION
    chfPressureTable_.set
    (
        new scalarFieldFieldInterpolateTable
        (
            chfPressureValues_,
            chfData_,
            chfPressureMethod_
        )
    );
    #elif OPENFOAMESI
    chfPressureTable_.reset
    (
        new scalarFieldFieldInterpolateTable
        (
            chfPressureValues_,
            chfData_,
            chfPressureMethod_
        )
    );
    #endif

    Info << "Initializing GroeneveldTable CHF model..." << endl;
}


scalar chfModel::compute
(
    const scalar p,
    const scalar T0,
    const scalar Q,
    const scalar G,
    const scalar DH,
    const scalar qWall,
    const label  faceI
) const
{
    switch (model_)
    {
        case Model::EPRI:
        {
            if (mag(qWall) <= SMALL)
            {
                return ROOTVGREAT;
            }

            const scalar p1(0.5328);
            const scalar p2(0.1212);
            const scalar p3(1.6151);
            const scalar p4(1.4066);
            const scalar p5(-0.3040);
            const scalar p6(0.4843);
            const scalar p7(-0.3285);
            const scalar p8(-2.0749);

            // G in Mlb/(h·ft²)
            const scalar GMlb =
                G / (1e6 * LbToKg_/HrToSec_/pow(FeetToMeters_, 2.0));

            const scalar reducedPressure = p/22064000.0;
            const scalar Fa = pow(GMlb, 0.1);
            const scalar Fc = 1.183*Fa;

            const scalar A =
                Fa * p1 * pow(reducedPressure, p2)
              * pow(GMlb, p5 + p7*reducedPressure);

            const scalar C =
                Fc * FAP_[faceI] * p3 * pow(reducedPressure, p4)
              * pow(GMlb, p6 + p8*reducedPressure);

            // q'' in MBtu/(h·ft²)
            const scalar qMBtu =
                qWall / (1e6 * BtuToJ_/HrToSec_/pow(FeetToMeters_, 2.0));

            return
                (
                    (A - Qin_)
                  /
                    (
                        C
                      + (Q - Qin_) / qMBtu
                    )
                )
              * (1e6 * BtuToJ_/HrToSec_/pow(FeetToMeters_, 2.0));
        }

        case Model::GroeneveldTable:
        {
            if (!chfPressureTable_.valid())
            {
                FatalErrorInFunction
                    << "GroeneveldTable CHF model selected but table was not "
                    << "initialised. Call initializeTable() first."
                    << abort(FatalError);
            }

            FieldField<Field, scalar> pressureData
            (
                (*chfPressureTable_)(p/1000.0)
            );

            scalarFieldInterpolateTable massFluxTable
            (
                chfMassFluxValues_,
                pressureData,
                chfMassFluxMethod_
            );

            scalarField qualityData(massFluxTable(G));

            scalarInterpolateTable qualityTable
            (
                chfQualityValues_,
                qualityData,
                chfQualityMethod_,
                interpolateTableBase::outOfBoundsMethod::FIXED
            );

            scalar qCHF = qualityTable(Q) * 1000.0;

            if (chfPitchToDiameter_ != 0.0)
            {
                const scalar K2 =
                    (2.0*chfPitchToDiameter_ - 1.5)
                  * exp(-cbrt(max(Q, -0.5))/2.0);

                qCHF *= K2;
            }
            else
            {
                const scalar K1 = max(0.6, sqrt(0.008/DH));
                qCHF *= K1;
            }

            return qCHF;
        }

        case Model::Zuber:
        {
            const scalar Tsat = IF97::Tsat97(p);

            return
                0.131 * IF97::latentHeatVaporization_p(p)
              * pow(IF97::rhovap_p(p), 0.5)
              * pow
                (
                    IF97::sigma97(Tsat)
                  * g_
                  * (IF97::rholiq_p(p) - IF97::rhovap_p(p)),
                    0.25
                );
        }

        case Model::ModifiedZuber:
        {
            const scalar Tsat = IF97::Tsat97(p);

            const scalar fsubc =
                1.0
              + 0.1
              * pow(IF97::rholiq_p(p)/IF97::rhovap_p(p), 0.75)
              * IF97::cpliq_p(p)
              * (Tsat - T0)
              / IF97::latentHeatVaporization_p(p);

            return
                0.131 * fsubc
              * IF97::latentHeatVaporization_p(p)
              * pow(IF97::rhovap_p(p), 0.5)
              * pow
                (
                    IF97::sigma97(Tsat)
                  * g_
                  * (IF97::rholiq_p(p) - IF97::rhovap_p(p)),
                    0.25
                );
        }

        case Model::GE:
        {
            if (mag(qWall) <= SMALL)
            {
                return ROOTVGREAT;
            }

            // G in lb/(h·ft²)
            const scalar GLb =
                G / (LbToKg_/HrToSec_/pow(FeetToMeters_, 2.0));

            const scalar unitConv = BtuToJ_/HrToSec_/pow(FeetToMeters_, 2.0);

            if (GLb < 500000.0)
            {
                return (840000.0 - 1000000.0*Q) * unitConv;
            }
            else
            {
                if (GLb >= 750000.0)
                {
                    WarningInFunction
                        << "Mass flux is out of range for GE CHF correlation."
                        << endl;
                }

                return (800000.0 - 1000000.0*Q) * unitConv;
            }
        }
    }

    FatalErrorInFunction
        << "Unhandled CHF model."
        << abort(FatalError);

    return ROOTVGREAT;
}

} // End namespace Foam

// ************************************************************************* //
