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

#include "rheology.H"
#include "zeroGradientFvPatchField.H"
#include "fvc.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rheology, 0);
    defineRunTimeSelectionTable(rheology, dictionary);
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rheology::rheology
(
    const fvMesh& mesh,    
    const materials& mat,
    const dictionary& rheologyDict
)
:     
    mesh_(mesh),
    mat_(mat),
    rheologyDict_(rheologyDict),
    epsTh_
    (
        IOobject
        (
            "epsTh",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("epsTh", dimless, symmTensor::zero)
    ),
    epsEl_
    (
        IOobject
        (
            "epsEl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("epsEl", dimless, symmTensor::zero)
    ),
    sigmaHyd_
    (
        IOobject
        (
            "sigmaHyd",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("sigmaHyd", dimPressure, 0.0)
    ),
    gradSigmaHyd_
    (
        IOobject
        (
            "gradSigmaHyd",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        fvc::grad(sigmaHyd_, "grad(sigmaHyd)")
    ),
    sigmaDev_
    (
        IOobject
        (
            "sigmaDev",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("sigmaDev", dimPressure, symmTensor::zero)
    ),
    additionalStrain_
    (
        IOobject
        (
            "additionalStrain",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("additionalStrain", dimless, symmTensor::zero),
        calculatedFvPatchField<scalar>::typeName
    ),
    thermalExpansion_(true),
    planeStress_(false),
    modifiedPlaneStrain_(false),
    solvePressureEqn_(false),
    pressureSmoothingScaleFactor_(1.0),
    residual_(0.0)
{                 
#ifdef OPENFOAMFOUNDATION
    if(!rheologyDict_.isNull())
#elif OPENFOAMESI        
    if(!rheologyDict_.isNullDict())
#endif
    {
        thermalExpansion_ = 
        rheologyDict_.lookupOrDefault<bool>("thermalExpansion", true);
    
        if(thermalExpansion_)
        {
            Info << nl << "     Thermal expansion: on" << nl;
        }
        else
        {
            Info <<  nl << "     Thermal expansion: off" << nl;
        }

        planeStress_ = 
        rheologyDict_.lookupOrDefault<bool>("planeStress", false);

        if (planeStress_)
        {
            Info<< "     Plane stress: on" << nl;
            Info<< "     Modified plane strain: off\n" << endl;
        }  
        else
        {
            //- Activate planeStain only if planeStress is off
            modifiedPlaneStrain_ = 
            rheologyDict_.lookupOrDefault<bool>("modifiedPlaneStrain", false);
            if(modifiedPlaneStrain_)
            {
                Info<< "     Plane stress: off" << nl;
                Info<< "     Modified plane strain: on\n" << endl;          
            }
            else
            {
                Info<< "     Plane stress: off" << nl;
                Info<< "     Modified plane strain: off\n" << endl;          
            }
        }

        solvePressureEqn_ = 
        rheologyDict_.lookupOrDefault<bool>("solvePressureEqn", false);

        pressureSmoothingScaleFactor_ = 
        rheologyDict_.lookupOrDefault<scalar>
        (
            "pressureSmoothingScaleFactor", 1.0
        );
    }
}


Foam::autoPtr<Foam::rheology>
Foam::rheology::New
(
    const fvMesh& mesh,
    const materials& mat,
    const dictionary& solverDict
)
{
    // Initialize type for mechanicsSolver class
    word type;

    dictionary rheologyDict
    (
        solverDict.subOrEmptyDict("rheologyOptions")
    );

    solverDict.lookup("rheology") >> type;
    
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn("rheology::New(const fvMesh&, const materials&, const dictionary&)")
            << "Unknown rheology type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<rheology>
    (
        cstrIter()
        (
            mesh, mat, rheologyDict
        )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rheology::~rheology()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::rheology::solvePressureEqn(volScalarField& sigmaHydExplicit)
{
    const volScalarField& mu(mat_.mu());
    const volScalarField& lambda(mat_.lambda());

    // Store previous iteration to allow relaxation
    sigmaHyd_.storePrevIter();

    // Lookup the momentum equation inverse diagonal field
    const volScalarField* ADPtr = NULL;

    if (mesh_.foundObject<volScalarField>("DEqnA"))
    {
        ADPtr = &mesh_.lookupObject<volScalarField>("DEqnA");
    }
    else if (mesh_.foundObject<volScalarField>("DDEqnA"))
    {
        ADPtr = &mesh_.lookupObject<volScalarField>("DDEqnA");
    }
    else
    {
        FatalErrorIn
        (
            "void " + type() + "updateSigmaHyd(...)\n"
        )   << "Cannot find the DEqnA or DDEqnA field: this should be "
            << "stored in the mechanicsSubSolver" << abort(FatalError);
    }

    const volScalarField& AD = *ADPtr;

    // Pressure diffusivity field
    const surfaceScalarField rDAf
    (
        "rDAf",
        pressureSmoothingScaleFactor_*fvc::interpolate
        (
            (2*mu+lambda)/AD, "interpolate(rDAf)"
        )
    );
    const dimensionedScalar one("one", dimless, 1.0);

    // Solve pressure laplacian
    // Note: the fvm and fvc laplacian terms cancel at convergence and the
    // laplacian - div(grad) term produce a smoothing/diffusion to quell
    // oscillations
    fvScalarMatrix sigmaHydEqn
    (
        fvm::Sp(one, sigmaHyd_)
      - fvm::laplacian(rDAf, sigmaHyd_, "laplacian(rDA,sigmaHyd)")
     ==
        sigmaHydExplicit
      - fvc::div(rDAf*fvc::interpolate(gradSigmaHyd_) & mesh_.Sf())
    );
        
    //- Compute residual
    residual_ = sigmaHydEqn.solve().max().initialResidual();    

    // Relax the pressure field
    sigmaHyd_.relax();

    gradSigmaHyd_ = fvc::grad(sigmaHyd_, "grad(sigmaHyd)");
}

 
bool Foam::rheology::converged()
{
    bool converged(false);

    const dictionary& stressControl = 
    mesh_.solutionDict().subDict("stressAnalysis");

    scalar convergenceTolerance(stressControl.lookupOrDefault("sigmaHyd",1e-6));

    //- Consider converge if both initialResidual and relResidual are 
    //- below threshold or if one of the two is 10 times below threshold.
    converged = ( 
                    (residual() < convergenceTolerance)
                );

    return converged;
}

// ************************************************************************* //
