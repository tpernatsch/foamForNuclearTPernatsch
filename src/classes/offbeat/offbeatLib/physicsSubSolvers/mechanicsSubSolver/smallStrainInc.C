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

#include "smallStrainInc.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smallStrainInc, 0);
    addToRunTimeSelectionTable
    (
        mechanicsSubSolver,
        smallStrainInc,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smallStrainInc::smallStrainInc
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
    )
{
    // Store old time fields
    D_.oldTime();
    gradD_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smallStrainInc::~smallStrainInc()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::smallStrainInc::correct()
{
    const dictionary& stressControl =
    mesh_.solutionDict().subDict("stressAnalysis");

    // Maximum inner iterations
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

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

    // Inner loop convergence variables
    int nInnerIter(0);
    bool convergedInner(false);

    do
    {
        // Store DD field for under-relaxation and residual calculation
        DD_.storePrevIter();
        D_.storePrevIter();

         // Update explicit component of the stress
        sigmaExp_ =
        mechanicsSubSolver::correctedSigmaExp(sigma_, twoMuLambda, gradDD_);

        // Update the divergence of the explicit component of the stress
        divSigmaExp_ =
        mechanicsSubSolver::correctedDivSigmaExp
        (sigmaExp_, sigma_, twoMuLambda, twoMuLambda_f, DD_, gradDD_);

        // Define incremental displacement equation
        fvVectorMatrix DDEqn
        (
            fvm::d2dt2(mat_.rho(), DD_)
          + fvc::d2dt2(mat_.rho().oldTime(), D_.oldTime())
        ==
            fvmLaplacianD(twoMuLambda_f, DD_, "laplacian(DD,D)")
          + divSigmaExp_
        );

        // Manipulate matrix based on BC rule (necessary for implicit contact)
        DDEqn.boundaryManipulate(DD_.boundaryFieldRef());

        forAll(refPairs, i)
        {
            DDEqn.setReference(refPairs[i].first(), refPairs[i].second());
        }

        // Relax matrix
        DDEqn.relax();

        // Solve equation for DD
        residual_ = DDEqn.solve().initialResidual();

        // Store initialResidual when in first inner iteration
        if(nInnerIter==0)
        {
            // Rescale residuals:
            initialResidual_ = residual_;
        }

        // Relax field
        DD_.relax();

        // Update the total displacement
        D_ = D_.oldTime() + DD_;

        // Update strain
        gradDD_ = mechanicsSubSolver::correctedGradD(twoMuLambda, DD_, gradDD_, sigmaExp_);
        gradD_ = gradD_.oldTime() + gradDD_;

        epsilon_ = symm(gradD_);

        // Save matrix coefficient
        const volScalarField DDEqnA("DDEqnA", DDEqn.A());

        // Update stress
        rheo_.correct(sigma(), epsilon(), DD());

        // Update explicit component of the stress
        sigmaExp_ = 
        mechanicsSubSolver::correctedSigmaExp(sigma_, twoMuLambda, gradDD_);
        sigma_ = symm(sigmaExp_ + twoMuLambda*gradDD_);

        // Calculate a different residual based on the relative change of DD
        scalar denom = max(gMax(mag(DD_.primitiveField())), SMALL);

        relResidual_ =
        gMax
        (
            mag(DD_.primitiveField() - DD_.prevIter().internalField())
        )/denom;

        Info << "relResidualD is " << relResidual_ << endl;

        nInnerIter++;

        convergedInner = converged(DD(), residual_);

    } while
    (
        not(convergedInner)
        && nInnerIter < nCorr
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
