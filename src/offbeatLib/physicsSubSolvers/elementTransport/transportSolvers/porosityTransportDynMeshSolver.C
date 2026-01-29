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

#include "porosityTransportDynMeshSolver.H"
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
#include "fvMesh.H"




// * * * * * * * * * * * * * Static Data Member  * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(porosityTransportDynMeshSolver, 0);
    addToRunTimeSelectionTable
    (
        porosityTransportDynMeshSolver, 
        porosityTransportDynMeshSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::porosityTransportDynMeshSolver::setDiffCoeff()
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
        } 
    }

    d_.correctBoundaryConditions();
}

void Foam::porosityTransportDynMeshSolver::readHolePatches()
{
    // Get list of patch names 
    wordList names(0);
    (porosityDict_.lookup("patches")) >> names;

    holePatchIDs_.setSize(names.size());
    
    forAll(names, i)
    {
        label patchID(dynMesh_.boundary().findPatchID(names[i]));

        if( patchID == -1)
        {
            FatalErrorInFunction()
            << "The name \"" << names[i] << "\" given in the \"patches\" subDict "
            << "in the porosityOptions dict" << nl
            << "does not correspond to any registered patch name." << nl << nl
            << "Possible patch are:" << nl << nl;
            forAll(dynMesh_.boundary(), j)
            {
                FatalErrorInFunction() << dynMesh_.boundaryMesh()[j].name() << endl;
            }
            FatalErrorInFunction() << abort(FatalError);
        }

        // Save addressing
        holePatchIDs_[i] = patchID;
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porosityTransportDynMeshSolver::porosityTransportDynMeshSolver
(
    fvMesh& mesh,
    const materials& mat,
    const dictionary& elementTransportDict,
    const word& solverName
)
:
    transportSolver(mesh, mat, elementTransportDict, solverName),
    dynMesh_(dynamic_cast<dynamicFvMesh&>(mesh)),
    T_(dynMesh_.lookupObject<volScalarField>("T")),
    gradT_(dynMesh_.lookupObject<volVectorField>("gradT")),
    poreVelocity_(dynMesh_.lookupObjectRef<volVectorField>("poreVelocity")),
    porosity_(dynMesh_.lookupObjectRef<volScalarField>("porosity")),
    d_
    (
        IOobject
        (
            "poreDiffusionCoeff",
            dynMesh_.time().timeName(),
            dynMesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        dynMesh_,
        dimensionedScalar("", dimArea/dimTime, 0.0),
        "zeroGradient"
    ),
    porosityDict_(elementTransportDict.subOrEmptyDict("porosityOptions")),
    diffusiveTerm_(porosityDict_.lookupOrDefault<bool>("diffusiveTerm", true)),
    boundedAdvectionTerm_(porosityDict_.lookupOrDefault<bool>("boundedAdvectionTerm", false)),
    updateMesh_(porosityDict_.lookupOrDefault<bool>("updateMesh", true)),
    holePatchIDs_(),
    maxCoAllowed_(GREAT),
    vCoalParameters_(),
    moveHole_(porosityDict_.lookupOrDefault<bool>("moveHole", true)),
    holePosition_(porosityDict_.lookupOrDefault<vector>("holePosition", vector(0,0,0))),
    holeMoved_(false),
    holeMoving_(false),
    useRiCorrection_(porosityDict_.lookupOrDefault<bool>("useRiCorrection", true)),
    useDirectionPV_(porosityDict_.lookupOrDefault<bool>("useDirectionPV", false))
{
    if(updateMesh_)
    {
#ifdef OPENFOAMFOUNDATION        
        // Read central hole patch IDs
        readHolePatches();

#elif OPENFOAMESI
        FatalErrorInFunction()
        << "The update mesh capabilities of porosityTransportDynMeshSolver "
        << "are not yet available for OpenFOAM ESI version." << nl
        << "Deactivate updateMesh." << abort(FatalError);
#endif
    }

    // Set diffusion coefficient 
    setDiffCoeff();

    // Set coalescence parameters
    if ( updateMesh_ and porosityDict_.found("vCoalescenceParameters")  )
    {
        dictionary parDict
        (
            porosityDict_.subDict("vCoalescenceParameters")
        );
        vCoalParameters_.setSize(parDict.size());
        vCoalParameters_[0] = readScalar(parDict.lookup("a"));
        vCoalParameters_[1] = readScalar(parDict.lookup("b"));
        vCoalParameters_[2] = readScalar(parDict.lookup("c"));
        vCoalParameters_[3] = readScalar(parDict.lookup("d"));
        vCoalParameters_[4] = readScalar(parDict.lookup("k"));
    }
    else if (updateMesh_ and ! porosityDict_.found("vCoalescenceParameters"))
    {
        FatalErrorInFunction()
        << "vCoalescenceParameters subDict need to be provided in "
        << "porosityOptions dictionary if updateMesh is activated." 
        << abort(FatalError);
        
        // vCoalParameters_.setSize(5);
        // vCoalParameters_[0] = 0.0;
        // vCoalParameters_[1] = 0.0;
        // vCoalParameters_[2] = 0.0;
        // vCoalParameters_[3] = 0.0;
        // vCoalParameters_[4] = 0.0;
    }

    // Read time stepping criterion
    const dictionary& controlDict(dynMesh_.time().controlDict());

    bool adjustTime = false;
    if (controlDict.found("adjustableTimeStep"))
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustableTimeStep", false);
    }
    else
    {
        adjustTime = controlDict.lookupOrDefault<bool>(
            "adjustTimeStep", false);
    }

    if(adjustTime)
    {
        controlDict.lookup("maxCourantNoPorosity") >> maxCoAllowed_;
    }

    // Store old time fields
    porosity_.oldTime();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
Foam::vector Foam::porosityTransportDynMeshSolver::computeHoleCenter()
{
    // Compute current center of the hole
    scalar nFaces = 0;
    vector sumCenters(0,0,0);
    forAll(holePatchIDs_, i)
    {
        const label patchID(holePatchIDs_[i]);
        const fvPatch& patch(dynMesh_.boundary()[patchID]);

        forAll(patch.Cf(), j)
        {
            sumCenters += patch.Cf()[j];
            nFaces += 1;
        }
    }

    return sumCenters/nFaces;
}

void Foam::porosityTransportDynMeshSolver::moveInitialHole()
{
    // Get hole current center
    const vector currentCenter = computeHoleCenter();

    // --- Compute new hole position ---

    // // Search maximum fuel temperature cell coordinates
    // scalar maxT(-1);
    // label cellMin(-1);

    // const labelListList& matAddrList(mat_.matAddrList());

    // forAll(mat_.materialsList(), i)
    // {
    //     if(isA<fuelMaterial>(mat_.materialsList()[i]))
    //     {
    //         const fuelMaterial& fuelMat = 
    //             refCast<const fuelMaterial>(mat_.materialsList()[i]);

    //         const dictionary matDict = fuelMat.materialModelDict();

    //         const labelList& addr(matAddrList[i]);

    //         forAll(addr, j)
    //         {
    //             const label cellI = addr[j];

    //             if(T_[cellI] > maxT)
    //             {
    //                 maxT = T_[cellI];
    //                 cellMin = cellI;
    //             }    
    //         }
    //     }
    // }

    // // Get coordinates of max T
    // vector newPosition = dynMesh_.C()[cellMin];

    // Read hole position
    vector newPosition = holePosition_;

    // --- Compute new hole position (end) ---

    // Get ref to pointDisplacement field (used to set BC for motionSolver)
    pointVectorField& pointDisp
    (
        dynMesh_.lookupObjectRef<pointVectorField>("pointDisplacement")
    );

    // Get ref to cellDisplacement (used for old-time disp at moving patch)
    const volVectorField& cellDisp
    (
        dynMesh_.lookupObjectRef<volVectorField>("cellDisplacement")
    );

    forAll(holePatchIDs_, i)
    {
        // Get current patch ID
        const label holePatchID(holePatchIDs_[i]);
        
        // Get ref to fvPatch
        const fvPatch& patch(dynMesh_.boundary()[holePatchID]);

        // Setup interpolator for this patch
        primitivePatchInterpolation patchInterpolator
        (
            dynMesh_.boundaryMesh()[holePatchID]
        );

        // Compute displacement at patch faces
        vectorField dispAtFaces(patch.size(), newPosition - currentCenter); 

        // Get old-time displacement (face value)
        const vectorField& oldDispFaces(cellDisp.boundaryField()[holePatchID]);

        // Interpolate to obtain displacement at patch points
        vectorField dispAtPoints = 
        patchInterpolator.faceToPointInterpolate
        (
            oldDispFaces + dispAtFaces
        );

        // Assign values to cellDisplacement (actually pointDisplacement) boundaryField
        pointDisp.boundaryFieldRef()[holePatchID] == dispAtPoints;
    }

    // Set bool to avoid repeating this function
    holeMoved_ = true;

    dynMesh_.update();

    // Print some info
    Info << "Current center of the hole: " << currentCenter << endl;
    Info << "New center of the hole: " << newPosition << nl << endl;
}

const Foam::scalarField Foam::porosityTransportDynMeshSolver::coalescenceVelocity
(
    const label patchID
)
{
    // Asymptotic porosity
    const scalar pLim(0.75);
    // const scalar pLim(1.00);

    // Unpack v coalescence parameters
    const scalar aaa(vCoalParameters_[0]);
    const scalar bbb(vCoalParameters_[1]);
    const scalar ccc(vCoalParameters_[2]);
    const scalar ddd(vCoalParameters_[3]);
    const scalar kkk(vCoalParameters_[4]);

    // Get ref to fvPatch
    const fvPatch& patch(dynMesh_.boundary()[patchID]);

    scalarField p(patch.size(), 0.0);

    // Get porosity on this patch
    p = (porosity_.boundaryField()[patchID]);

    // Compute gradient of pore velocity normal to patch
    const vectorField& pvPatch(poreVelocity_.boundaryField()[patchID]);
    const vectorField pvPatchCells(poreVelocity_.internalField(), patch.faceCells());
    scalarField gradPVp
    (
        ((pvPatch - pvPatchCells)*patch.deltaCoeffs()) & patch.nf()

    );

    // Get inner radii for this patch 
    scalarField r_i(patch.size(), 1.0);
    if( useRiCorrection_ )
    {
        // Get current hole center
        // const vector hc = computeHoleCenter();
        const vector hc(0,0,0);
        forAll(r_i, i)
        {
            r_i[i] = 
                pow
                ( pow(patch.Cf()[i].x()-hc.x(), 2.0) 
                + pow(patch.Cf()[i].y()-hc.y(), 2.0) 
                , 0.5 );
        }
    } 

    // With limiter
    scalarField coalVel =  
    max
    (
        // // gradV, r_i, 2nd deg poly, asymptotic
        // kkk *(mag(gradPVp) * pow(r_i, ddd) * (aaa * pow(min(p, pLim),2.0) + bbb * min(p, pLim) + ccc))/(pLim-min(p, 0.99*pLim))
        
        // gradV, r_i, exponential, asymptotic
        kkk * (mag(gradPVp) * pow(r_i, ddd) * (aaa * (exp(bbb * min(p, pLim)) - 1.0)))/(pLim-min(p, 0.99*pLim))

        ,

        0.0
    );

    return coalVel;
}

void Foam::porosityTransportDynMeshSolver::updateMesh()
{
    if ( updateMesh_ )
    {
        // Initialize central void position
        // TODO: to be checked with restart
        if ( moveHole_ and not holeMoved_ )
        {
            Info << nl << "Initializing position of central hole ..." << endl;

            // Set to true to prevent fvc::makeRelative
            holeMoving_ = true;

            // Move initial hole in the mesh
            moveInitialHole();

            dynMesh_.moving(false);
        }
        else
        {
            // Set to false to allow fvc::makeRelative
            holeMoving_ = false;

            // Get ref to pointDisplacement field (used to set BC for motionSolver)
            pointVectorField& pointDisp
            (
                dynMesh_.lookupObjectRef<pointVectorField>("pointDisplacement")
            );

            // Get ref to cellDisplacement (used for old-time disp at moving patch)
            const volVectorField& cellDisp
            (
                dynMesh_.lookupObjectRef<volVectorField>("cellDisplacement")
            );

            forAll(holePatchIDs_, i)
            {
                // Get current patch ID
                const label holePatchID(holePatchIDs_[i]);
                
                // Get ref to fvPatch
                const fvPatch& patch(dynMesh_.boundary()[holePatchID]);

                // Setup interpolator for this patch
                primitivePatchInterpolation patchInterpolator
                (
                    dynMesh_.boundaryMesh()[holePatchID]
                );

                // Compute coalescence velocity for this patch
                const scalarField coalVel = coalescenceVelocity(holePatchID);

                // Compute displacement at patch faces
                vectorField dispAtFaces = vectorField(patch.size(), vector::zero);
                if( useDirectionPV_ )
                {
                    vectorField dir(patch.size(), vector::zero);
                    forAll(dir, i_dir)
                    {
                        const vector& pv = poreVelocity_[patch.faceCells()[i_dir]];
                        dir[i_dir] = pv/max(mag(pv),VSMALL);
                    }
                    // dispAtFaces = - coalVel * dynMesh_.time().deltaTValue() * dir;
                    dispAtFaces = - (coalVel * dynMesh_.time().deltaTValue() * dir & patch.nf()) * patch.nf();
                }
                else
                {
                    dispAtFaces = - coalVel * dynMesh_.time().deltaTValue() * patch.nf();
                }

                // Get old-time displacement (face value)
                const vectorField& oldDispFaces(cellDisp.boundaryField()[holePatchID]);

                // Interpolate to obtain displacement at patch points
                vectorField dispAtPoints = 
                patchInterpolator.faceToPointInterpolate
                (
                    oldDispFaces + dispAtFaces
                );

                // When dispAtFaces is very small, the interpolation at points
                // can do mistakes and net result is a tiny displacement in 
                // opposite direction wrt dispAtFaces, that cumulated over time
                // creates a closing (instead of opening) of the hole  
                if( max(mag(dispAtFaces)) > 1e-10 )
                {
                    // Assign values to cellDisplacement (actually pointDisplacement) boundaryField
                    pointDisp.boundaryFieldRef()[holePatchID] == dispAtPoints;
                }
            }

            // Update the dynamic mesh
            dynMesh_.update();
            
            // Compute total porosity balance - post mesh update
            scalar poreVolNew(gSum(porosity_*dynMesh_.V().field()));

            Info << "Porosity info:" << nl
                << tab << "Total pore volume : " << poreVolNew << nl
                << endl;
        }
    }
}

void Foam::porosityTransportDynMeshSolver::correct()
{
    const dictionary& stressControl = 
        dynMesh_.solutionDict().subDict("stressAnalysis");
    
    const int nCorr = stressControl.lookupOrDefault<int>("nCorrectors", 1);

    int nInnerIter = 0;

    bool convergedInner(false);

    do
    {
        // Store prev iter values for rel residual
        porosity_.storePrevIter();

        // Flux for advective term
        surfaceScalarField flux = fvc::interpolate( poreVelocity_ ) &  dynMesh_.Sf();

        // Include mesh velocity in advective flux
        // if( not holeMoving_ )
        // {
            fvc::makeRelative(flux, poreVelocity_);
        // }

        // Equation for pore transport (only advection contribution)
        fvScalarMatrix poreEqn
        (
              fvm::ddt(porosity_)
            + fvm::div(flux, porosity_, "div(flux,P)")
        );

        // Multiplies the porosity in div term by (1-p)
        if ( boundedAdvectionTerm_ )
        {
            poreEqn -= fvc::div(flux, pow(porosity_,2.0));
        }

        // Add diffusive term to the equation
        if( diffusiveTerm_ )
        {
            poreEqn -= fvm::laplacian(d_, porosity_);
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

    Info << endl;
}

Foam::scalar Foam::porosityTransportDynMeshSolver::nextDeltaT()
{
    // Flux for Courant number
    const surfaceScalarField phi = 
        fvc::interpolate((poreVelocity_)*porosity_) & dynMesh_.Sf();

    // Define tmp field of Courant numbers
    tmp<volScalarField> tCourantNo
    (
        volScalarField::New
        (
            "CourantNo",
            dynMesh_,
            dimensionedScalar(dimless, 0),
            "zeroGradient"
        )
    );

    // Compute Courant number (according to OpenFOAM method)
    tCourantNo->ref() = (0.5*dynMesh_.time().deltaT())*fvc::surfaceSum(mag(phi))()()/dynMesh_.V();
    tCourantNo->correctBoundaryConditions();

    // Compute maximum Courant number 
    const scalar maxCo(gMax(tCourantNo->field()));

    scalar currentDeltaT = dynMesh_.time().deltaT().value();

    // Compute maximum allowed deltaT multiplier 
    scalar maxDtMultiplier = maxCoAllowed_/max(maxCo, SMALL);

    // Compute maximum allowed next deltaT
    scalar nextDeltaT( currentDeltaT * maxDtMultiplier );

    Info << "Maximum deltaT calculated by porosityTransportDynMeshSolver: " 
         << dynMesh_.time().timeToUserTime(nextDeltaT)
         << " with a maximum Courant No of " << maxCo 
         << endl;

    return min(nextDeltaT, 1e10);
}

// ************************************************************************* //
