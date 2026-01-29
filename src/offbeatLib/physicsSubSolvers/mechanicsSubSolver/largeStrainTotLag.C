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

#include "largeStrainTotLag.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"
#include "logVolFields.H"
#include "accelerationSchemes.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(largeStrainTotLag, 0);
    addToRunTimeSelectionTable
    (
        mechanicsSubSolver, 
        largeStrainTotLag, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::largeStrainTotLag::largeStrainTotLag
(
    fvMesh& mesh,
    const materials& mat,
    rheology& rheo,
    const dictionary& mechanicsDict
)
:    
    mechanicsSubSolver(mesh, mat, rheo, mechanicsDict),
    strainTensorName_("EulerAlmansi"),
    F_
    (
        IOobject
        (
            "F",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("I", dimless, I)
    ),
    Finv_
    (
        IOobject
        (
            "Finv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(F_)
    ),
    relF_
    (
        IOobject
        (
            "relF",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        F_ & inv(F_.oldTime())
    ),
    J_
    (
        IOobject
        (
            "J",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(F_)
    ),
    storedPrevIter_(false)
{
    // Store old time fields
    D_.oldTime();
            
    if(
#ifdef OPENFOAMFOUNDATION
        !mechanicsDict_.isNull()
#elif OPENFOAMESI        
        !mechanicsDict_.isNullDict()
#endif
    )
    {
        strainTensorName_ = 
        mechanicsDict_.lookupOrDefault<word>("strainTensor", "EulerAlmansi");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::largeStrainTotLag::~largeStrainTotLag()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::largeStrainTotLag::correct() 
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

    // Intermediate fields
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
        // Previous-previous iteration for inertial relaxation
        DispPrevPrevIter = D_.prevIter();

        // Store fields for under-relaxation and residual calculation
        D_.storePrevIter(); 

        // Update explicit component of the stress
        sigmaExp_ = 
        mechanicsSubSolver::correctedSigmaExp(
            (J_*(Finv_ & sigma_))(), twoMuLambda, gradD_);

        // Update the divergence of the explicit component of the stress
        divSigmaExp_ = 
        mechanicsSubSolver::correctedDivSigmaExp
        (sigmaExp_, (J_*(Finv_ & sigma_))(), twoMuLambda, twoMuLambda_f, D_, gradD_);
        
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
        
        // Un-couple implicit contact boundaries (if present)
        coupleBoundaries(D_, false);

        // Relax field
        D_.relax();      

        // Apply acceleration scheme
        if (accelScheme.valid())
        {
            accelScheme->accelerate();
        }

        // Update strain and stress
        gradD_ = 
        mechanicsSubSolver::correctedGradD(twoMuLambda, D_, gradD_, sigmaExp_);

        // Total deformation gradient
        F_ = I + gradD_.T();
        
        // Inverse of the deformation gradient
        Finv_ = inv(F_);

        // Relative deformation gradient
        relF_ = F_ & inv(F_.oldTime());

        // Jacobian of the deformation gradient
        J_ = det(F_);

        // Calculate epsilon 
        if(strainTensorName_ ==  "linear")
        {
            epsilon_ = symm(gradD_);
        }
        else if(strainTensorName_ ==  "EulerAlmansi")
        {
            epsilon_ = 0.5*(I - symm(Finv_.T() & Finv_));
        }
        else if(strainTensorName_ ==  "Hencky")
        {
            epsilon_ = 0.5*log(symm(F_.T() & F_));
        }
        else if(strainTensorName_ ==  "GreenLagrange")
        {
            epsilon_ = 0.5*(symm(F_.T() & F_) - I);
        }
        else
        {
            // Error for strain tensor
            FatalErrorInFunction() << nl  
            << "Possible strain tensors are:" << nl  
            << tab << "- linear" << nl
            << tab << "- EulerAlmansi" << nl
            << tab << "- Hencky" << nl
            << tab << "- GreenLagrange" << nl
            << exit(FatalError);
        }

        // Save matrix coefficient
        const volScalarField DEqnA("DEqnA", DEqn.A());

        // Correct constitutive material laws
        rheo_.correct(sigma(), epsilon(), D());

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

        convergedInner = mechanicsSubSolver::converged(D(), residual_, relResidual_, absErr_);
        
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


