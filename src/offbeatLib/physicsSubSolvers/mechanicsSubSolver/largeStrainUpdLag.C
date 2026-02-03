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

#include "largeStrainUpdLag.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "addToRunTimeSelectionTable.H"
#include "processorPolyPatch.H"
#include "twoDPointCorrector.H"
#include "symmetryPolyPatch.H"
#include "logVolFields.H"
#include "accelerationSchemes.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(largeStrainUpdLag, 0);
    addToRunTimeSelectionTable
    (
        mechanicsSubSolver, 
        largeStrainUpdLag, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::largeStrainUpdLag::largeStrainUpdLag
(
    fvMesh& mesh,
    const materials& mat,
    rheology& rheo,
    const dictionary& mechanicsDict
)
:     
    mechanicsSubSolver(mesh, mat, rheo, mechanicsDict),
    strainTensorName_("EulerAlmansi"),
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
        I + gradDD_.T()
    ),
    relFinv_
    (
        IOobject
        (
            "relFinv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        inv(relF_)
    ),
    relJ_
    (
        IOobject
        (
            "relJ",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        det(relF_)
    ),
    referenceMesh_
    (        
       IOobject
       (
           "referenceMesh",
           mesh.time().constant(),
           mesh
       ),
       pointField
       (
         pointIOField
         (
            IOobject
            (
                "points",
                mesh.time().constant()/polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
         )
       ),
       faceList(mesh.faces()),
       cellList(mesh.cells())
    ),
    initialPoints_
    (
        IOobject
        (
            "initialPoints",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh.points()
    ),
    storedPrevIter_(false)
{
    IOobject D_header( "D", mesh_.time().timeName(), mesh_, IOobject::NO_READ);
    if (D_header.typeHeaderOk<volVectorField>()) {
        Warning << "Using incremental solver with initial displacement field \"D\". "
                << "Are you re-starting a simulation?" << endl;
    }

    // Store old time fields
    D_.oldTime();
    DD_.oldTime();
    F_.oldTime();
    J_.oldTime();
    gradD_.oldTime();
            
    if
    (
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

Foam::largeStrainUpdLag::~largeStrainUpdLag()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::largeStrainUpdLag::correct()
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
        sigmaExp_ = 
        mechanicsSubSolver::correctedSigmaExp(
            (relJ_*(relFinv_ & sigma_))(), twoMuLambda, gradDD_);

        // Update the divergence of the explicit component of the stress
        divSigmaExp_ = 
        mechanicsSubSolver::correctedDivSigmaExp
        (sigmaExp_, (relJ_*(relFinv_ & sigma_))(), twoMuLambda, twoMuLambda_f, DD_, gradDD_);        

        // Couple implicit contact boundaries (if present)
        coupleBoundaries(DD_, true);

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
        gradDD_ = mechanicsSubSolver::correctedGradD(twoMuLambda, DD_, gradDD_, sigmaExp_);
        
        // Calculate gradient of total displacement
        gradD_ = gradD_.oldTime() + gradDD_;

        // Relative deformation gradient
        relF_ = I + gradDD_.T();

        // Inverse relative deformation gradient
        relFinv_ = inv(relF_);

        // Total deformation gradient
        F_ = relF_ & F_.oldTime();

        // Inverse of total deformation gradient
        Finv_ = inv(F_);

        // Relative Jacobian (Jacobian of relative deformation gradient)
        relJ_ = det(relF_);

        // Jacobian of deformation gradient
        J_ = relJ_*J_.oldTime();

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
            << tab << "- GreenLagrange"  << nl
            << exit(FatalError);
        }

        // Save matrix coefficient
        const volScalarField DDEqnA("DDEqnA", DDEqn.A());
        
        // Correct constitutive material laws
        rheo_.correct(sigma(), epsilon(), DD());
        
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


