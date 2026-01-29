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
    const dictionary& thermalDict
)
:
    thermalSubSolver(mesh, mat, thermalDict),
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
        if(mesh_.foundObject<volScalarField>("Q"))
        {
            Q_ = &mesh_.lookupObject<volScalarField>("Q");
        }
    }
    
    // Maximum inner iterations
    const int nCorrectors = thermalSubSolver::nCorrectors(T_.name());

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
        );

        if(Q_ != nullptr)
        {
            // Optional heat source.
            // Since the equation is built without an explicit RHS, the source term
            // is applied by subtracting it from the LHS:  LHS = Q  ⇔  LHS - Q = 0.
            TEqn -= *Q_;
        }
        
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
        }        

        // Calculate a different residual based on the relative change of T
        scalar denom = 
        gMax(mag(T_.primitiveField() - T_.oldTime().primitiveField()));

        absErr_ = 
        gMax(mag(T_.primitiveField() - T_.prevIter().internalField()));

        if (denom < SMALL)
        {
            denom = max(gMax(mag(T_.primitiveField())), SMALL);
        }
        
        relResidual_ = absErr_/denom;

        if(denom < 1e-5)
        {
            relResidual_ = VSMALL;
        }

        Info << "relResidualT is " << relResidual_ << endl;
        Info << "absErrT is " << absErr_ << endl;

        nInnerIter++;

        convergedInner = converged();

    } while 
    (
        not(convergedInner)  
        && nInnerIter < nCorrectors 
    );

    Info << endl;
}

// ************************************************************************* //
