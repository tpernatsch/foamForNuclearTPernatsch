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
    DispPrev(D_.size(), Foam::vector::zero),
    DispPrevPrev(D_.size(), Foam::vector::zero),
    inertialRelaxation_(mechanicsDict.lookupOrDefault<bool>("inertialRelaxation", false)),
    inertialRelaxationDeltaT_(mechanicsDict.lookupOrDefault<scalar>("inertialRelaxationDeltaT", 1e-7))
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

    do
    {
        // Store fields for under-relaxation and residual calculation
        DispPrevPrev = DispPrev;
        D_.storePrevIter(); 
        DispPrev = D_.prevIter();


        // Update explicit component of the stress
        sigmaExp_ = 
        mechanicsSubSolver::correctedSigmaExp(sigma_, twoMuLambda, gradD_);

        // Update the divergence of the explicit component of the stress
        divSigmaExp_ = 
        mechanicsSubSolver::correctedDivSigmaExp
        (sigmaExp_, sigma_, twoMuLambda, twoMuLambda_f, D_, gradD_);

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
            DEqn.setReference(refPairs[i].first(), refPairs[i].second());
        }


        //-Possible relaxation technique derived from ddt2
        if( inertialRelaxation_ )
        { 
            scalar deltaT  = inertialRelaxationDeltaT_;
            scalar deltaT0 = inertialRelaxationDeltaT_;
         
            scalar coefft   = (deltaT + deltaT0)/(2*deltaT);
            scalar coefft00 = (deltaT + deltaT0)/(2*deltaT0);
            scalar coefft0  = coefft + coefft00;
         
            scalar rDeltaT2 = 4.0/sqr(deltaT + deltaT0);
         
            {
                DEqn.diag() += (coefft*rDeltaT2)*mesh_.V()*mat_.rho().primitiveField();
        
                DEqn.source() += rDeltaT2*mesh_.V()*mat_.rho().primitiveField()*
                (
                    coefft0*DispPrev
                  - coefft00*DispPrevPrev
                );
            }
        }

        
        // Relax matrix
        DEqn.relax();
        
        // Solve equation for displacement
        residual_ = DEqn.solve().initialResidual();  

        // Store initialResidual when in first inner iteration
        if(nInnerIter==0)
        {
            initialResidual_ = residual_; 
        }

        // Relax field
        D_.relax();     

        // Update strain 
        gradD_ = 
        mechanicsSubSolver::correctedGradD(twoMuLambda, D_, gradD_, sigmaExp_);
        epsilon_ = symm(gradD_);

        // Save matrix coefficient
        const volScalarField DEqnA("DEqnA", DEqn.A());     
        
        // Update stress
        rheo_.correct(sigma(), epsilon(), D());   

        // Update explicit component of the stress
        sigmaExp_ = 
        mechanicsSubSolver::correctedSigmaExp(sigma_, twoMuLambda, gradD_);
        sigma_ = symm(sigmaExp_ + twoMuLambda*gradD_);

        // Calculate a different residual based on the relative change of D
        scalar denom = gMax
        (
            mag(D_.primitiveField() - DOld.primitiveField())
        );

        if (denom < SMALL)
        {
            denom = max(gMax(mag(D_.primitiveField())), SMALL);
        }
        
        relResidual_ =
        gMax
        (
            mag(D_.primitiveField() - D_.prevIter().internalField())
        )/denom;

        Info << "relResidualD is " << relResidual_ << endl;

        nInnerIter++;

        convergedInner = converged(D(), residual_);

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


