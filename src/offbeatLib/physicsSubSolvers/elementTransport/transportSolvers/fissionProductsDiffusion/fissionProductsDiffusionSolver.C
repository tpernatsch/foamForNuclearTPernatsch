/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fissionProductsDiffusionSolver.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"
#include "UO2.H"
#include "wordRe.H"
#include "uniformDimensionedFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fissionProductsDiffusionSolver, 0);
    addToRunTimeSelectionTable
    (
        transportSolver,
        fissionProductsDiffusionSolver,
        dictionary
    );

    static const scalar eJPerFission = 3.2e-11;   // J per fission event
    static const scalar Na           = 6.0221408e23; // Avogadro constant
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fissionProductsDiffusionSolver::updateCoefficients()
{
    scalarField Bu(mesh_.nCells(), scalar(0));
    if (Bu_) Bu = scalarField(Bu_->primitiveField());

    const scalarField  f  = Q_.primitiveField() / eJPerFission;

    const labelListList& matAddrList = mat_.matAddrList();

    forAll(FpNames_, k)
    {
        // Yield Y [atoms/fission] — Q_=0 in non-fissioning zones, so s=0 there
        scalarField Y(mesh_.nCells(), scalar(0));
        yieldModels_[k]->updateYield(Y, Bu);
        s_[k]().primitiveFieldRef() = Y * f / Na;

        // Optional external source [mol/m³/s] injected by the application
        // (e.g. compact–TRISO coupling). The field must be named
        // "extSrc_N_<species>" and registered on this mesh.
        const word extSrcName("extSrc_N_" + FpNames_[k]);
        if (mesh_.foundObject<volScalarField>(extSrcName))
        {
            s_[k]().primitiveFieldRef() +=
                mesh_.lookupObject<volScalarField>(extSrcName).internalField();
        }

        // Per-zone diffusion coefficient
        forAll(mat_.materialsList(), i)
        {
            if (diffModels_[k][i].valid())
            {
                diffModels_[k][i]->updateCoef
                (
                    diff_[k](), T_, matAddrList[i], FpNames_[k]
                );
            }
        }

        diff_[k]().correctBoundaryConditions();
    }
}


void Foam::fissionProductsDiffusionSolver::outputInfo()
{
    const scalarField volume(mesh_.V());
    const scalar dT = mesh_.time().deltaT().value();

    forAll(FpNames_, i)
    {
        Fpleft_[i]    = gSum(volume * Fp_[i]());
        totalProd_[i] = totalProdPrev_[i]
                      + gSum(volume * s_[i]()) * dT;

        const scalar inventoryOld =
            gSum(volume * Fp_[i]().oldTime().internalField());

        const scalar deltaReleased =
            inventoryOld - Fpleft_[i] + gSum(volume * s_[i]()) * dT;

        releasedInc_[i] = deltaReleased / max(dT, SMALL);
        released_[i]    = releasedPrev_[i] + deltaReleased;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fissionProductsDiffusionSolver::fissionProductsDiffusionSolver
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:
    transportSolver(mesh, mat, elementTransportDict, solverName),
    T_(mesh_.lookupObject<volScalarField>("T")),
    phi_
    (
        mesh_.foundObject<volScalarField>("fastFluence")
        ? &mesh_.lookupObject<volScalarField>("fastFluence") : nullptr
    ),
    Bu_
    (
        mesh_.foundObject<volScalarField>("Bu")
        ? &mesh_.lookupObject<volScalarField>("Bu") : nullptr
    ),
    Q_(mesh_.lookupObject<volScalarField>("Q")),
    releasedInc_(),
    releasedPrev_(),
    released_(),
    timeIndex_(-1),
    totalInitial_(),
    totalProd_(),
    totalProdPrev_(),
    Fpleft_(),
    enrichment_(0),
    FpNames_()
{
    const dictionary FpDict =
        elementTransportDict.subOrEmptyDict("FpDiffusionOptions");

    // Species names — word list: fissionProducts (Cs Ag Kr ...);
    FpDict.lookup("fissionProducts") >> FpNames_;

    // Read enrichment from UO2 material (needed by enrichmentDependent yield model)
    forAll(mat_.materialsList(), i)
    {
        if (isA<UO2>(mat_.materialsList()[i]))
        {
            enrichment_ =
                refCast<const UO2>(mat_.materialsList()[i]).enrichment();
            break;
        }
    }

    // Allocate per-species fields
    Fp_.resize(FpNames_.size());
    diff_.resize(FpNames_.size());
    s_.resize(FpNames_.size());

    forAll(FpNames_, i)
    {
        Fp_[i].set
        (
            new volScalarField
            (
                IOobject
                (
                    "N_" + FpNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimMoles/dimVolume, 0.0),
                "fixedValue"
            )
        );

        Fp_[i]().storeOldTime();


        diff_[i].set
        (
            new volScalarField
            (
                IOobject
                (
                    "diffN_" + FpNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimArea/dimTime, 0.0),
                "calculated"
            )
        );

        s_[i].set
        (
            new volScalarField
            (
                IOobject
                (
                    "srcN_" + FpNames_[i],
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar("", dimMoles/dimVolume/dimTime, 0.0),
                "calculated"
            )
        );
    }

    const label nSpec = FpNames_.size();
    const label nMat  = mat_.materialsList().size();

    releasedInc_.resize(nSpec, 0);
    releasedPrev_.resize(nSpec, 0);
    released_.resize(nSpec, 0);
    totalInitial_.resize(nSpec, 0);
    totalProd_.resize(nSpec, 0);
    totalProdPrev_.resize(nSpec, 0);
    Fpleft_.resize(nSpec, 0);

    // Build per-species per-zone diffusion coefficient models.
    // diffusion {} is material-keyed; zone patterns (wordRe) are the keys.
    // A per-species subdict inside a zone dict overrides the zone-level type.
    const dictionary& diffDict = FpDict.subDict("diffusion");
    const wordList diffPatterns = diffDict.toc();

    diffModels_.resize(nSpec);
    forAll(FpNames_, k)
    {
        diffModels_[k].resize(nMat);
        forAll(mesh_.cellZones(), i)
        {
            const word& zoneName = mesh_.cellZones()[i].name();
            forAll(diffPatterns, p)
            {
                wordRe wre(diffPatterns[p], wordRe::DETECT);
                if (wre.match(zoneName))
                {
                    const dictionary& zoneDict = diffDict.subDict(diffPatterns[p]);
                    // Per-species subdict takes priority over zone-level "type"
                    const dictionary& modelDict =
                        zoneDict.found(FpNames_[k])
                        ? zoneDict.subDict(FpNames_[k])
                        : zoneDict;
                    const word modelType =
                        modelDict.lookupOrDefault<word>("type", word::null);
                    diffModels_[k][i] =
                        diffCoefModel::New(mesh_, modelDict, modelType);
                    break;
                }
            }
        }
    }

    // Build per-species yield models.
    // yield {} is material-keyed; only fissile zones need entries.
    // A per-species subdict overrides the zone-level "type".
    // The solver injects "fpName" into the dict so smart-default models
    // (e.g. type UO2) can dispatch internally per species.
    const dictionary& yieldDictRoot = FpDict.subOrEmptyDict("yield");
    const wordList yieldPatterns = yieldDictRoot.toc();

    yieldModels_.resize(nSpec);
    forAll(FpNames_, k)
    {
        bool found = false;
        forAll(yieldPatterns, p)
        {
            if (found) break;
            forAll(mesh_.cellZones(), i)
            {
                wordRe wre(yieldPatterns[p], wordRe::DETECT);
                if (wre.match(mesh_.cellZones()[i].name()))
                {
                    const dictionary& zoneDict = yieldDictRoot.subDict(yieldPatterns[p]);
                    const dictionary& modelDict =
                        zoneDict.found(FpNames_[k])
                        ? zoneDict.subDict(FpNames_[k])
                        : zoneDict;
                    dictionary modelDictWithFp(modelDict);
                    modelDictWithFp.set("fpName", FpNames_[k]);
                    yieldModels_[k] = fpYieldModel::New(modelDictWithFp, enrichment_);
                    found = true;
                    break;
                }
            }
        }
        if (!found)
            FatalErrorInFunction
                << "No yield model found for species " << FpNames_[k]
                << " in yield {} dict."
                << abort(FatalError);
    }

    // Initial inventory
    const scalarField volume(mesh_.V());
    forAll(FpNames_, k)
    {
        totalInitial_[k] = gSum(volume * Fp_[k]());
        Fp_[k]().correctBoundaryConditions();
    }

    updateCoefficients();

    // Register uniformDimensionedScalarFields per species for function-object access
    forAll(FpNames_, i)
    {
        // Release fraction (dimensionless)
        {
            uniformDimensionedScalarField* ptr = new uniformDimensionedScalarField
            (
                IOobject
                (
                    "fpRelease_" + FpNames_[i],
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar(dimless, scalar(0))
            );
            ptr->store(ptr);
        }
        // Cumulative released amount [mol] — for external coupling (e.g. global FGR)
        {
            uniformDimensionedScalarField* ptr = new uniformDimensionedScalarField
            (
                IOobject
                (
                    "releasedN_" + FpNames_[i],
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dimensionedScalar(dimMoles, scalar(0))
            );
            ptr->store(ptr);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fissionProductsDiffusionSolver::correctFp
(
    volScalarField& Fp,
    const volScalarField& diffFp,
    const volScalarField& s
)
{
    label nInnerIter = 0;
    const label nCorr = transportSolver::nCorrectors(Fp.name());

    bool convergedInner = false;
    do
    {
        Fp.storePrevIter();

        fvScalarMatrix FpEqn
        (
              fvm::ddt(Fp)
            - fvm::laplacian(diffFp, Fp, "laplacianFp")
         ==
            s
        );

        FpEqn.relax();
        residual_ = FpEqn.solve().max().initialResidual();

        if (nInnerIter == 0)
        {
            initialResidual_ = residual_;
        }

        Fp.relax();

        scalar denom =
            gMax(mag(Fp.primitiveField() - Fp.oldTime().primitiveField()));

        const scalar num =
            gMax(mag(Fp.primitiveField() - Fp.prevIter().internalField()));

        if (denom < SMALL)
        {
            denom = max(gMax(mag(Fp.primitiveField())), SMALL);
        }

        relResidual_ = num/denom;

        Info << "relResidual " << Fp.name() << " " << relResidual_ << endl;

        convergedInner = converged(Fp.name());

    } while (!convergedInner && ++nInnerIter < nCorr);
}


void Foam::fissionProductsDiffusionSolver::correct()
{
    updateCoefficients();

    if (mesh_.time().timeIndex() != timeIndex_)
    {
        releasedPrev_ = released_;
        totalProdPrev_ = totalProd_;
        timeIndex_    = mesh_.time().timeIndex();
    }

    scalarList residualFp(FpNames_.size());
    forAll(FpNames_, i)
    {
        correctFp(Fp_[i](), diff_[i](), s_[i]());
        residualFp[i] = initialResidual_;
    }

    initialResidual_ = max(residualFp);

    outputInfo();

    Info << "Fission products transport:" << endl;
    forAll(FpNames_, i)
    {
        const scalar releaseFraction = min
        (
            max
            (
                released_[i] / max(VSMALL, totalProd_[i] + totalInitial_[i]),
                scalar(0)
            ),
            scalar(1)
        );

        // Keep registered scalars current so function objects and the app can read them
        mesh_.lookupObjectRef<uniformDimensionedScalarField>
        (
            "fpRelease_" + FpNames_[i]
        ).value() = releaseFraction;

        mesh_.lookupObjectRef<uniformDimensionedScalarField>
        (
            "releasedN_" + FpNames_[i]
        ).value() = released_[i];

        Info << tab << "Initial "   << FpNames_[i] << " amount: "    << totalInitial_[i] << nl
             << tab << "Production of " << FpNames_[i] << ": "       << totalProd_[i]   << nl
             << tab << "Released "  << FpNames_[i] << " amount: "    << released_[i]    << nl
             << tab << "Remained "  << FpNames_[i] << " amount: "    << Fpleft_[i]      << nl
             << tab << "Fraction of " << FpNames_[i] << " release: " << releaseFraction << nl;
    }
    Info << endl;
}

// ************************************************************************* //
