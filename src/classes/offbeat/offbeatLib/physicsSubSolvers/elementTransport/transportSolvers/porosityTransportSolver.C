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

#include "porosityTransportSolver.H"
#include "fvm.H"
#include "fvc.H"
#include "typeInfo.H"
#include "addToRunTimeSelectionTable.H"
#include "mechanicsSubSolver.H"
#include "layerAdditionRemovalPolyTopoChanger.H"
#include "polyTopoChange.H"
#include "primitivePatchInterpolation.H"
#include "PrimitivePatchInterpolation.H"
#include "wedgePolyPatch.H"



// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porosityTransportSolver, 0);
    addToRunTimeSelectionTable
    (
        porosityTransportSolver, 
        porosityTransportSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porosityTransportSolver::setDiffCoeff()
{
    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        if(isA<fuelMaterial>(mat_.materialsList()[i]))
        {
            const fuelMaterial& fuelMat = 
                refCast<const fuelMaterial>(mat_.materialsList()[i]);

            const dictionary matDict = fuelMat.materialModelDict();

            const scalar d = 
                matDict.lookupOrDefault<scalar>("poreDiffusionCoeff", 1e-12);

            const labelList& addr(matAddrList[i]);

            forAll(addr, j)
            {
                const label cellI = addr[j];

                d_[cellI] = d;
            }
            // Test: diffusion only for cells adjacent to hole patches
            // forAll(holePatchIDs_, hpid)
            // {
            //     const polyPatch& patch = mesh_.boundaryMesh()[hpid];
            //     forAll(patch.faceCells(), fc)
            //     {
            //         d_[patch.faceCells()[fc]] = d;
            //     }
            // }
        }
    }

    d_.correctBoundaryConditions();
}


void Foam::porosityTransportSolver::updateRestructuringFields()
{
    // Update central void field
    centralVoid_ *= 0.0;
    forAll(porosity_, i)
    {
        if(porosity_[i] >= voidThreshold_)
        {
            centralVoid_[i] = 1;
        }
    }

    // Update restructured region field
    restructuredRegion_ *= 0.0;

    // Compute porosity gradient
    volScalarField gradP = mag(fvc::grad(porosity_));

    forAll(gradP, i)
    {
        const scalar porosityChange = 
            mag(porosity_[i] - asFabricatedP_[i])/max(asFabricatedP_[i],VSMALL); 
        if
        ( 
            porosityChange > 0.01
            or
            gradP[i] > 1
        )
        {
            restructuredRegion_[i] = 1.0;
        }
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityTransportSolver::porosityTransportSolver
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:
    transportSolver(mesh, mat, elementTransportDict, solverName),
    T_(mesh_.lookupObject<volScalarField>("T")),
    gradT_(mesh_.lookupObject<volVectorField>("gradT")),
    poreVelocity_(mesh_.lookupObject<volVectorField>("poreVelocity")),
    porosity_(mesh_.lookupObjectRef<volScalarField>("porosity")),
    d_
    (
        IOobject
        (
            "poreDiffusionCoeff",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimArea/dimTime, 0.0),
        "zeroGradient"
    ),
    centralVoid_
    (
        IOobject
        (
            "centralVoid",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        "zeroGradient"
    ),
    restructuredRegion_
    (
        IOobject
        (
            "restructuredRegion",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        "zeroGradient"
    ),
    asFabricatedP_
    (
        IOobject
        (
            "asFabricatedP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("", dimless, 0.0),
        "zeroGradient"
    ),
    porosityDict_(elementTransportDict.subOrEmptyDict("porosityOptions")),
    diffusiveTerm_(porosityDict_.lookupOrDefault<bool>("diffusiveTerm", false)),
    boundedAdvectionTerm_(porosityDict_.lookupOrDefault<bool>("boundedAdvectionTerm", false)),
    updateMesh_(porosityDict_.lookupOrDefault<bool>("updateMesh", false)),
    voidThreshold_(porosityDict_.lookupOrDefault<scalar>("voidThreshold", 1.0)),
    maxCoAllowed_(GREAT)
{
    // Set diffusion coefficient 
    setDiffCoeff();

    // Read time stepping criterion
    const dictionary& cDict(mesh_.time().controlDict());

    bool adjustTime(cDict.lookup("adjustableTimeStep"));

    if(adjustTime)
    {
        cDict.lookup("maxCourantNoPorosity") >> maxCoAllowed_;
    }

    // Store as fabricated (i.e. initial) porosity 
    // might cause problems for restrart
    asFabricatedP_ = porosity_;
    asFabricatedP_.correctBoundaryConditions();

    // Store old time fields
    porosity_.oldTime();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::porosityTransportSolver::correct()
{

    const dictionary& stressControl = 
        mesh_.solutionDict().subDict("stressAnalysis");
    
    const int nCorr = 
        stressControl.lookupOrDefault<int>("nCorrectors", 1);

    int nInnerIter = 0;

    bool convergedInner(false);

    do
    {
        // Store prev iter values for rel residual
        porosity_.storePrevIter();

        // Flux for advective term
        surfaceScalarField flux = fvc::interpolate( poreVelocity_ ) &  mesh_.Sf();

        // Equation for pore transport (only advection contribution)
        fvScalarMatrix poreEqn
        (
              fvm::ddt(porosity_)
            + fvm::div(flux, porosity_, "div(flux,P)")
        );

        // Add diffusive term to the equation
        if( diffusiveTerm_ )
        {
            poreEqn -= fvm::laplacian(d_, porosity_);
        }

        // Multiplies the porosity in div term by (1-p)
        if ( boundedAdvectionTerm_ )
        {
            poreEqn -= fvc::div(flux, pow(porosity_,2.0));
        }
        
        // Relax equation    
        poreEqn.relax();
        
        // Compute residual
        residual_ = poreEqn.solve().max().initialResidual();

        // Store initial residual if in first inner iteration
        if(nInnerIter==0)
        {
            initialResidual_ = residual_; 
        }
        
        // Relax fields
        porosity_.relax();

        // Calculate a different residual based on the relative change
        scalar denom = 
        gMax(mag(porosity_.primitiveField() - porosity_.oldTime().primitiveField()));

        scalar num = 
        gMax(mag(porosity_.primitiveField() - porosity_.prevIter().internalField()));

        if (denom < SMALL)
        {
            denom = max(gMax(mag(porosity_.primitiveField())), SMALL);
        }
        
        relResidual_ = num/denom;

        if(denom < 1e-5)
        {
            relResidual_ = VSMALL;
        }

        Info << "relResidualPorosity " << relResidual_ << endl;

        convergedInner = converged(porosity_.name());

    } while 
    (
        not(convergedInner)  
        && ++nInnerIter < nCorr 
    );

    // Update centralVoid and restructuredRegion fields
    updateRestructuringFields();

    Info << endl;
}

Foam::scalar Foam::porosityTransportSolver::nextDeltaT()
{
    // Flux for Courant number
    const surfaceScalarField phi = 
        fvc::interpolate((poreVelocity_)*porosity_) & mesh_.Sf();

    // Define tmp field of Courant numbers
    tmp<volScalarField> tCourantNo
    (
        volScalarField::New
        (
            "CourantNo",
            mesh_,
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        )
    );

    // Compute Courant number (according to OpenFOAM method)
    tCourantNo->ref() = (0.5*mesh_.time().deltaT())*fvc::surfaceSum(mag(phi))()()/mesh_.V();
    tCourantNo->correctBoundaryConditions();

    // Compute maximum Courant number 
    const scalar maxCo(gMax(tCourantNo->field()));

    scalar currentDeltaT = mesh_.time().deltaT().value();

    // Compute maximum allowed deltaT multiplier 
    scalar maxDtMultiplier = maxCoAllowed_/max(maxCo, SMALL);

    // Compute maximum allowed next deltaT
    scalar nextDeltaT( currentDeltaT * maxDtMultiplier );

    Info << "Maximum deltaT calculated by porosityTransportSolver: " 
         << mesh_.time().timeToUserTime(nextDeltaT)
         << " with a maximum Courant No of " << maxCo 
         << endl;

    return min(nextDeltaT, 1e10);
}

// ************************************************************************* //
