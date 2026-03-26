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

#include "smallStrainIncUpdated.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"
#include "accelerationSchemes.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smallStrainIncUpdated, 0);
    addToRunTimeSelectionTable
    (
        mechanicsSubSolver,
        smallStrainIncUpdated,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smallStrainIncUpdated::smallStrainIncUpdated
(
    fvMesh& mesh,
    const materials& mat,
    rheology& rheo,
    const dictionary& mechanicsDict
)
:
    mechanicsSubSolver(mesh, mat, rheo, mechanicsDict),
    DD_
    (
        IOobject
        (
            "DD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("DD", dimLength, vector::zero),
        calculatedFvPatchField<scalar>::typeName
    ),
    gradDD_
    (
        IOobject
        (
            "gradDD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("gradDD", dimless, tensor::zero),
        calculatedFvPatchField<scalar>::typeName
    ),
    DSigma_
    (
        IOobject
        (
            "DSigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("DSigma", dimPressure, symmTensor::zero)
    ),
    DSigmaExp_
    (
        IOobject
        (
            "DSigmaExp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("DSigmaExp", dimPressure, tensor::zero)
    ),
    divDSigmaExp_
    (
        IOobject
        (
            "divDSigmaExp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("divDSigmaExp", dimPressure/dimLength, vector::zero)
    ),
    updateMesh_(mechanicsDict.lookupOrDefault<Switch>("updateMesh", true)),
    storedPrevIter_(false)
{

    IOobject D_header( "D", mesh_.time().timeName(), mesh_, IOobject::NO_READ);
    if (D_header.typeHeaderOk<volVectorField>()) {
        Warning << "Using incremental solver with initial displacement field \"D\". "
                << "Are you re-starting a simulation?" << endl;
    }

    if(updateMesh_)
    {
        referenceMesh_.reset
        (
            new fvMesh
            (
                IOobject
                (
                    "referenceMesh",
                    mesh_.time().constant(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                pointField
                (
                    pointIOField
                    (
                        IOobject
                        (
                            "points",
                            mesh_.time().constant()/polyMesh::meshSubDir,
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        )
                    )
                ),
                faceList(mesh_.faces()),
                cellList(mesh_.cells())
            )
        );

        // Write reference mesh
        referenceMesh_->write();
    }

    // Store old time fields
    DD_.oldTime();
    D_.oldTime();
    gradD_.oldTime();
    sigma_.oldTime();
    sigmaExp_.oldTime();
    divSigmaExp_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smallStrainIncUpdated::~smallStrainIncUpdated()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::smallStrainIncUpdated::updateIncrementalFields()
{
    gradDD_ = gradDD_*0;
    DD_ = DD_*0;
    DSigma_ = DSigma_*0;
}

void Foam::smallStrainIncUpdated::correct()
{
    const dictionary& stressControl =
    mesh_.solutionDict().subDict("stressAnalysis");

    // Maximum inner iterations
    const int nCorrectors = mechanicsSubSolver::nCorrectors(D_.name());

    // Acceleration scheme
    autoPtr<vectorAccelerationScheme> accelScheme;

    if (stressControl.found("accelerationScheme"))
    {
        accelScheme = vectorAccelerationScheme::New
        (
            1,
            stressControl.subDict("accelerationScheme")
        );
        accelScheme->addField(DD_);

        if (nCorrectors < accelScheme->minIter())
        {
            IOWarningInFunction(stressControl)
                << "Number of inner iterations for field \"" << D_.name() << "\" is too small. "
                << "Expected at least " << accelScheme->minIter() << " \"nCorrectors\" for this field. "
                << "Solution acceleration will have no effect. " << nl << endl;
        }
    }

    // Apply multi-material corrections
    if (multiMaterialCorrection_)
    {
        materialInterface_->correct();
    }

    List<Tuple2<label, vector> > refPairs;

    if (stressControl.found("referencePairs"))
    {
        refPairs = List<Tuple2<label, vector> >(stressControl.lookup("referencePairs"));
    }

    // Intermediate fields
    const volScalarField& mu(mat_.mu());
    const volScalarField& lambda(mat_.lambda());
    volScalarField twoMuLambda("2mu+lambda", 2*mu + lambda);
    surfaceScalarField twoMuLambda_f
        = mechanicsSubSolver::correctedTwoMuLambdaF(twoMuLambda);

    // Reference to old time
    // const volVectorField& DDOld(DD_.oldTime());

    // Inner loop convergence variables
    int nInnerIter(0);
    bool convergedInner(false);

    // Previous-previous iteration for inertial relaxation
    vectorField DispPrevPrevIter;
    if (!storedPrevIter_)
    {
        DD_.storePrevIter();
        storedPrevIter_ = true;
    }

    do
    {
        // Previous-previous iteration for inertial relaxation
        DispPrevPrevIter = DD_.prevIter();

        // Store DD field for under-relaxation and residual calculation
        DD_.storePrevIter();

        // Update explicit component of the stress
        DSigmaExp_ =
        mechanicsSubSolver::correctedSigmaExp(DSigma_, twoMuLambda, gradDD_);

        sigmaExp_ = sigma_.oldTime() + DSigmaExp_;

        // Couple implicit contact boundaries (if present)
        coupleBoundaries(DD_, true);

        // Define incremental displacement equation
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(mat_.rho(), DD_)
        ==
            fvmLaplacianD(twoMuLambda_f, DD_, "laplacian(DD,D)")
        );

        // Update the divergence of the explicit component of the stress
        divDSigmaExp_ =
        mechanicsSubSolver::correctedDivSigmaExp
        (DSigmaExp_, DSigma_, twoMuLambda, twoMuLambda_f, DD_, gradDD_);

        divSigmaExp_ = divSigmaExp_.oldTime() + divDSigmaExp_;

        DDEqn -= divDSigmaExp_;

        // Manipulate matrix based on BC rule (necessary for implicit contact)
        DDEqn.boundaryManipulate(DD_.boundaryFieldRef());

        forAll(refPairs, i)
        {
            DDEqn.setReference(refPairs[i].first(), refPairs[i].second(), true);
        }

        // Relax matrix
        mechanicsSubSolver::relax
        (
            DDEqn,
            DD_.prevIter(),
            DispPrevPrevIter
        );

        // Solve equation for DD
        residual_ = DDEqn.solve().initialResidual();

        // Un-couple implicit contact boundaries (if present)
        coupleBoundaries(DD_, false);

        // Relax field
        DD_.relax();

        // Apply acceleration scheme
        if (accelScheme.valid())
        {
            accelScheme->accelerate();
        }

        // Update the total displacement
        D_ = D_.oldTime() + DD_;

        // Calculate gradient of incremental displacement
        gradDD_ = mechanicsSubSolver::correctedGradD(twoMuLambda, DD_, gradDD_, DSigmaExp_);

        // Calculate gradient of total displacement
        gradD_ = gradD_.oldTime() + gradDD_;

        // Calculate strain
        epsilon_ = symm(gradD_);

        // Save matrix coefficient
        const volScalarField DDEqnA("DDEqnA", DDEqn.A());

        // Update stress
        volSymmTensorField sigma0 = sigma_.oldTime();
        rheo_.correct(sigma(), epsilon(), DD());
        sigmaEq_ = sqrt((3.0/2.0)*magSqr(dev(sigma_)));
        DSigma_ = sigma() - sigma0;
        
        // Calculate a different residual based on the relative change of D
        scalar denom = max(gMax(mag(DD_.primitiveField())), SMALL);

        absErr_ = gMax
        (
            mag(DD_.primitiveField() - DD_.prevIter().internalField())
        );
        
        relResidual_ = absErr_/denom;

        Info << "relResidualD is " << relResidual_ << endl;
        Info << "absErrD is " << absErr_ << endl;

        convergedInner = mechanicsSubSolver::converged(DD(), residual_, relResidual_, absErr_);
        
        // Store initialResidual when in first inner iteration
        if(nInnerIter==0)
        {
            initialResidual_ = residual_; 
            initialRelResidual_ = relResidual_; 
            initialAbsErr_ = absErr_; 
        }
        
        nInnerIter++;    

    } while
    (
        not(convergedInner)
        && nInnerIter < nCorrectors
    );

    // Transform stress and displacement in cyl. coord. (for output)
    correctCylindrical();
    // Transform stress and displacement in sph. coord. (for output)
    correctSpherical();

    // Print total force on each face
    forceSummary();

    Info<< "Domain average displacement: "
        << fvc::domainIntegrate(D_).value()/gSum(mesh_.V())
        << endl << endl;
}

// ************************************************************************* //
