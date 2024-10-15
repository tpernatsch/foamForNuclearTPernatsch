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

#include "solidConductionSolver.H"
#include "zeroGradientFvPatchField.H"
#include "calculatedFvPatchField.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solidConductionSolver, 0);
    addToRunTimeSelectionTable
    (
        thermalSubSolver, 
        solidConductionSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidConductionSolver::solidConductionSolver
(
    const fvMesh& mesh,
    materials& mat,
    const dictionary& thermalOptDict
)
:
    thermalSubSolver(mesh, mat, thermalOptDict),
    Q_(nullptr),
    hGap_
    (
        IOobject
        (
            "hGap",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("hGap", dimless, 0)
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::solidConductionSolver::correct()
{
    if(Q_ == nullptr)
    {
        Q_ = &mesh_.lookupObject<volScalarField>("Q");
    }
    
    const dictionary& stressControl = 
    mesh_.solutionDict().subDict("stressAnalysis");
    
    // Set maximum number of inner correction and residual tolerance
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);
    // scalar convergenceTolerance = readScalar(stressControl.lookup("T"));

    int nInnerIter = 0;

    bool convergedInner(false);

    do
    {
        T_.storePrevIter();
        
        //- Correct the thermoMech properties
        mat_.correctThermoMechProperties(); 
        
        volScalarField rhoCp("rhoCp", mat_.rho()*mat_.Cp());
        
        // Heat conduction equation
        fvScalarMatrix TEqn
        (
            fvm::ddt(rhoCp, T_)
           - fvm::laplacian(mat_.k(), T_)
            ==
            *Q_
        );
        
        // Relax equation    
        TEqn.relax();
        
        // Compute residual
        residual_ = TEqn.solve().max().initialResidual();

        // Compute gradT_
        computeGradT();

        if(calculateEnthalpy_)
        {
            calcEnthalpy();
        }
       
        // Store initial residual if in first inner iteration
        if(nInnerIter==0)
        {
            initialResidual_ = residual_; 
        }
        
        // Relax T field
        T_.relax();

        if (heatFluxSummary_)
        {
            calcHeatFlux(TEqn.flux());
            heatFluxForMappings(TEqn.flux());
        }        

        // Calculate a different residual based on the relative change of T
        scalar denom = 
        gMax(mag(T_.primitiveField() - T_.oldTime().primitiveField()));

        scalar num = 
        gMax(mag(T_.primitiveField() - T_.prevIter().internalField()));

        if (denom < SMALL)
        {
            denom = max(gMax(mag(T_.primitiveField())), SMALL);
        }
        
        relResidual_ = num/denom;

        if(denom < 1e-5)
        {
            relResidual_ = VSMALL;
        }

        Info << "relResidualT " << relResidual_ << endl;

        nInnerIter++;

        convergedInner = converged();

    } while 
    (
        not(convergedInner)  
        && nInnerIter < nCorr 
    );

    Info << endl;
}

// ************************************************************************* //
