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

#include "smallStrain.H"
#include "fvm.H"
#include "fvc.H"
#include "linear.H"
#include "harmonic.H"
#include "skewCorrectionVectors.H"
#include "zeroGradientFvPatchField.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"
#include "accelerationSchemes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smallStrain, 0);
    addToRunTimeSelectionTable
    (
        mechanicsSubSolver, 
        smallStrain, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smallStrain::smallStrain
(
    fvMesh& mesh,
    const materials& mat,
    rheology& rheo,
    const dictionary& mechanicsDict
)
:     
    mechanicsSubSolver(mesh, mat, rheo, mechanicsDict),
    storedPrevIter_(false)
{
    // Store old time fields
    D_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::smallStrain::~smallStrain()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::smallStrain::correct() 
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
        accelScheme->addField(D_);

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

    // Calculate twoMuLambda and twoMuLambdaf fields
    const volScalarField& mu(mat_.mu());
    const volScalarField& lambda(mat_.lambda());
    volScalarField twoMuLambda("2mu+lambda", 2*mu + lambda);
    surfaceScalarField twoMuLambda_f
        = mechanicsSubSolver::correctedTwoMuLambdaF(twoMuLambda);
    
    // Reference to old time
    const volVectorField& DOld(D_.oldTime());

    // Inner loop convergence variables
    int nInnerIter(0);
    bool convergedInner(false);

    // Previous-previous iteration for inertial relaxation
    vectorField DispPrevPrevIter;
    if (!storedPrevIter_)
    {
        D_.storePrevIter();
        storedPrevIter_ = true;
    }

    do
    {
        // Previous-previous iteration for equation relaxation
        DispPrevPrevIter = D_.prevIter();

        // Store fields for under-relaxation and residual calculation
        D_.storePrevIter();
        
        // Update explicit component of the stress
        sigmaExp_ = 
        mechanicsSubSolver::correctedSigmaExp(sigma_, twoMuLambda, gradD_);

        // Update the divergence of the explicit component of the stress
        divSigmaExp_ = 
        mechanicsSubSolver::correctedDivSigmaExp
        (sigmaExp_, sigma_, twoMuLambda, twoMuLambda_f, D_, gradD_);

        // Couple implicit contact boundaries (if present)
        coupleBoundaries(D_, true);

        // Define displacement equation
        fvVectorMatrix DEqn
        (
            fvm::d2dt2(mat_.rho(), D_)
          ==
            fvmLaplacianD(twoMuLambda_f, D_, "laplacian(DD,D)")
          + divSigmaExp_
        );
                            
        // Manipulate matrix based on BC rule (necessary for implicit contact)
        DEqn.boundaryManipulate(D_.boundaryFieldRef()); 

        forAll(refPairs, i)
        {
            DEqn.setReference(refPairs[i].first(), refPairs[i].second(), true);
        }

        // Relax matrix
        mechanicsSubSolver::relax
        (
            DEqn,
            D_.prevIter(),
            DispPrevPrevIter
        );
        
        // Solve equation for displacement
        residual_ = DEqn.solve().initialResidual();

        // Un-coupled implicit contact boundaries (if present)
        coupleBoundaries(D_, false);

        // Relax field
        D_.relax();

        // Apply acceleration scheme
        if (accelScheme.valid())
        {
            accelScheme->accelerate();
        }

        // Update strain 
        gradD_ = 
        mechanicsSubSolver::correctedGradD(twoMuLambda, D_, gradD_, sigmaExp_);
        epsilon_ = symm(gradD_);

        // Save matrix coefficient
        const volScalarField DEqnA("DEqnA", DEqn.A());     
        
        // Update stress
        rheo_.correct(sigma(), epsilon(), D());
        sigmaEq_ = sqrt((3.0/2.0)*magSqr(dev(sigma_)));

        // Calculate a different residual based on the relative change of D
        scalar denom = gMax
        (
            mag(D_.primitiveField() - DOld.primitiveField())
        );

        if (denom < SMALL)
        {
            denom = max(gMax(mag(D_.primitiveField())), SMALL);
        }

        absErr_ = gMax
        (
            mag(D_.primitiveField() - D_.prevIter().internalField())
        );
        
        relResidual_ = absErr_/denom;

        Info << "relResidualD is " << relResidual_ << endl;
        Info << "absErrD is " << absErr_ << endl;

        convergedInner = converged(D(), residual_, relResidual_, absErr_);
        
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


