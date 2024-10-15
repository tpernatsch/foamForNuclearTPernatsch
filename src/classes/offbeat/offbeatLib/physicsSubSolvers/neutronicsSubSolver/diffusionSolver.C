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

#include "diffusionSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(diffusionSolver, 0);
    addToRunTimeSelectionTable
    (
        neutronicsSubSolver, 
        diffusionSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusionSolver::diffusionSolver
(
    const fvMesh& mesh,
    const burnup& bu,
    const dictionary& neutronicsOptDict
)
:
    neutronicsSubSolver(mesh, bu, neutronicsOptDict),
    neutronFlux0_
    (
        IOobject
        (
            "neutronFlux0",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless/dimArea/dimTime, 1.0)
    ),
    invDiffLength_
    (
        IOobject
        (
            "invDiffLength",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless/dimArea, 0.0),
        "zeroGradient"
    ),
    SigmaA_
    (
        IOobject
        (
            "SigmaA",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless/dimLength, 0),
        "zeroGradient"
    ),
    neutronD_
    (
        IOobject
        (
            "neutronD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimLength, GREAT),
        "zeroGradient"
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::diffusionSolver::correct()
{
    const dictionary& neutronicsControl = 
    mesh_.solutionDict().subDict("stressAnalysis");
    
    const int nCorr = neutronicsControl.lookupOrDefault<int>("nCorrectors", 1);
    
    scalar convergenceTolerance(neutronicsControl.lookupOrDefault("neutronFlux0",1e-6));

    int nInnerIter = 0;

    // Reference to volScalarField internal fields
    scalarField& SigmaAI   = SigmaA_.ref();
    scalarField& neutronDI = neutronD_.ref();

    // Get SigmaA
    SigmaAI = burnup_.SigmaA();

    // Get neutronD
    neutronDI = burnup_.neutronD();
    neutronD_.correctBoundaryConditions();

    // Compute square inverse diffusion length
    invDiffLength_ = SigmaA_ / neutronD_;

    do
    {
        // Store prev iter for relResidual calculation
        neutronFlux0_.storePrevIter();

        // Thermal neutron diffusion equation
        fvScalarMatrix neutronDiffusionEqn
        (
              fvm::laplacian(neutronFlux0_)
            - fvm::Sp(invDiffLength_, neutronFlux0_)
        );
        
        // Relax equation    
        neutronDiffusionEqn.relax();
        
        // Compute residual
        residual_ = neutronDiffusionEqn.solve().max().initialResidual();    
       
        // Store initial residual if in first inner iteration
        if(nInnerIter==0)
        {
            initialResidual_ = residual_; 
        }
        
        // Relax field
        neutronFlux0_.relax();

        // Correct boundary values
        neutronFlux0_.correctBoundaryConditions();

        // Calculate a different residual based on the relative change of nFlux
        scalar denom = 
        gMax(mag(neutronFlux0_.primitiveField() - neutronFlux0_.oldTime().primitiveField()));

        scalar num = 
        gMax(mag(neutronFlux0_.primitiveField() - neutronFlux0_.prevIter().internalField()));

        if (denom < SMALL)
        {
            denom = max(gMax(mag(neutronFlux0_.primitiveField())), SMALL);
        }
        
        relResidual_ = num/denom;
        /*
        if(denom < 1e-5)
        {
            relResidual_ = VSMALL;
        }
        */

        Info << "relResidualneutronFlux0 " << relResidual_ << endl;

    } while (
                residual_ > convergenceTolerance 
                && 
                ++nInnerIter < nCorr
            );

    Info << endl;
}

 
bool Foam::diffusionSolver::converged()
{
    bool converged(false);

    const dictionary& stressControl = 
    mesh_.solutionDict().subDict("stressAnalysis");

    bool useRelRes(stressControl.lookupOrDefault("useRelResNeutronFlux0", true));
    
    scalar convergenceTolerance
    (
        stressControl.lookupOrDefault("neutronFlux0", 1e-6)
    );

    if ( useRelRes )
    {
        scalar relConvergenceTolerance
        (
            readScalar(stressControl.lookup("relneutronFlux0"))
        );

        // Consider converged if both initialResidual and relResidual are 
        // below threshold or if one of the two is 10 times below threshold.
        converged = 
        ( 
            initialResidual_ < convergenceTolerance 
            && 
            relResidual_ < relConvergenceTolerance 
        )
        or
        (
            initialResidual_ < convergenceTolerance/10
        )
        or
        (
            relResidual_ < relConvergenceTolerance/10 
        )
        ; 
    }
    else
    {
        // In this case only the residual computed by OpenFOAM is accounted for 
        // the convergence monitoring.
        converged = 
        ( 
            initialResidual_ < convergenceTolerance
        ); 
    }

    return converged;
}


// ************************************************************************* //
