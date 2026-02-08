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

#include "mechanicsSubSolver.H"
#include "leastSquaresGradf.H"
#include "fvm.H"
#include "fvc.H"
#include "calculatedFvPatchFields.H"
#include "zeroGradientFvPatchFields.H"
#include "calculatedTotalDisplacementFvPatchVectorField.H"
#include "tractionDisplacementFvPatchVectorField.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "wedgePolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "twoDPointCorrector.H"
#include "symmetryPolyPatch.H"
#include "surfaceInterpolation.H"
#include "harmonic.H"
#include "linear.H"
#include "snGradScheme.H"
#include "uncorrectedSnGrad.H"
#include "gaussLaplacianScheme.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mechanicsSubSolver, 0);
    defineRunTimeSelectionTable(mechanicsSubSolver, dictionary);    

    addToRunTimeSelectionTable
    (
        mechanicsSubSolver, 
        mechanicsSubSolver, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


void Foam::mechanicsSubSolver::relax
(
    fvVectorMatrix& fvm,
    const vectorField& psi0,
    const vectorField& psi00,
    const scalar alpha
)
{
    // Description:
    // This relaxation method is derived from fvMatrix::relax using the
    // idea of intertial (d2dt2) relaxation rather than pure matrix 
    // relaxation.
    // The standard fvMatrix relaxation method doesn't work well for the
    // displacement equation because it doesn't include the old-old
    // iteration value. The result is a sequence that almost never converges.
    // It was noticed, however, that the inertial method also fails to converge
    // if the old-old value is assumed zero. Consequently, we have updated
    // the original fvMatrix relaxation here to include the old-old value
    // while ensuring diagonal dominance. This approach seems to work
    // reasonably well for numerically ill-posed problems with free 
    // bodies.

    if (alpha <= 0)
    {
        return;
    }
    
    const volVectorField& psi = fvm.psi();

    // Matrix diagonal plus internal boundary coefficients
    scalarField D(fvm.diag());

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(D.size(), 0.0);
    fvm.sumMagOffDiag(sumOff);

    // Handle the boundary contributions to the diagonal
    forAll(psi.boundaryField(), patchi)
    {
        const fvPatchVectorField& ptf = psi.boundaryField()[patchi];

        if (ptf.size())
        {
            const labelUList& pa = fvm.lduAddr().patchAddr(patchi);
            vectorField& iCoeffs = fvm.internalCoeffs()[patchi];

            if (ptf.coupled())
            {
                const vectorField& pCoeffs = fvm.boundaryCoeffs()[patchi];

                forAll(pa, face)
                {
                    D[pa[face]] += mag(cmptMin(iCoeffs[face]));
                    sumOff[pa[face]] += cmptMax(cmptMag(pCoeffs[face]));
                }
            }
            else
            {
                forAll(pa, face)
                {
                    D[pa[face]] += mag(cmptMin(iCoeffs[face]));
                }
            }
        }
    }

    // Dominance ratio
    scalarField r = D / max(sumOff, SMALL);

    DebugInFunction
        << "Relaxing " << psi.name() << " by " << alpha << nl
        << "Minimum dominance ratio = " << gMin(r) << endl;

    r = min(alpha*r, 1.0);
    scalarField delta = (1 - r)/r * D;

    fvm.diag() += delta;
    fvm.source() += delta*(2*psi0 - psi00);
}


void Foam::mechanicsSubSolver::relax
(
    fvVectorMatrix& fvm,
    const vectorField& psi0,
    const vectorField& psi00
)
{
    const word& name = fvm.psi().name();
    const fvMesh& mesh = fvm.psi().mesh();

    if (mesh.relaxEquation(name))
    {
        relax(fvm, psi0, psi00, mesh.equationRelaxationFactor(name));
    }
}


void Foam::mechanicsSubSolver::updateMesh
(
    fvMesh& mesh,
    const volVectorField& movingD,
    const bool& keepFlatEmptyFaces
) 
{
    Info << "Updating mesh with " << movingD.name() << "..." << endl;

    //- Constant reference to class for mesh point interpolation
    const volPointInterpolation& meshPointInterpolation = 
        volPointInterpolation::New(mesh);

    //- Compute the displacement of the points (due to deltaD)
    pointVectorField pointsDisplacement = 
        meshPointInterpolation.interpolate(movingD);

    vectorField& pointsDisplacementI = pointsDisplacement.ref();

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<symmetryPolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            // const labelList& meshPoints =
            //     mesh().boundaryMesh()[patchI].meshPoints();

            // if
            // (
            //     returnReduce(mesh().boundaryMesh()[patchI].size(), sumOp<int>())
            // )
            // {
            //     continue;
            // }

            // const vector avgN =
            //     gAverage(mesh().boundaryMesh()[patchI].pointNormals());

            // const vector i(1, 0, 0);
            // const vector j(0, 1, 0);
            // const vector k(0, 0, 1);

            // if (mag(avgN & i) > 0.95)
            // {
            //     forAll(meshPoints, pI)
            //     {
            //         pointDDI[meshPoints[pI]].x() = 0;
            //     }
            // }
            // else if (mag(avgN & j) > 0.95)
            // {
            //     forAll(meshPoints, pI)
            //     {
            //         pointDDI[meshPoints[pI]].y() = 0;
            //     }
            // }
            // else if (mag(avgN & k) > 0.95)
            // {
            //     forAll(meshPoints, pI)
            //     {
            //         pointDDI[meshPoints[pI]].z() = 0;
            //     }
            // }
        }
        else if (
            isA<emptyPolyPatch>(mesh.boundaryMesh()[patchI])
            and keepFlatEmptyFaces)
        {
            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // if
            // (
            //     returnReduce(mesh().boundaryMesh()[patchI].size(), sumOp<int>())
            // )
            // {
            //     continue;
            // }

            const vector avgN =
                gAverage(mesh.boundaryMesh()[patchI].pointNormals());
                
            const vector k(0, 0, 1);

            if (mag(avgN & k) > 0.95)
            {
                forAll(meshPoints, pI)
                {
                    pointsDisplacementI[meshPoints[pI]].z() = 0;
                }
            }
        }
    }

    //- Calculate the new position of the displaced mesh points
    vectorField pointsDisplaced = 
        mesh.oldPoints()
    +   pointsDisplacement.internalField();

    twoDPointCorrector twoDCorrector(mesh);
    twoDCorrector.correctPoints(pointsDisplaced);

    //- Move the points of the mesh to the new positions
    mesh.movePoints(pointsDisplaced);

    mesh.moving(true);
#ifdef OPENFOAMFOUNDATION    
    mesh.setPhi().writeOpt() = IOobject::NO_WRITE;
#elif OPENFOAMESI    
    mesh.setPhi()->writeOpt() = IOobject::NO_WRITE;
#endif    

    //- Exclude solution in wedge direction if axisymmetric
    mechanicsSubSolver::checkWedges(mesh);

    Info << "Mesh updated.\n" << endl;
}


void Foam::mechanicsSubSolver::checkWedges(fvMesh& mesh)
{
/*********** Derived from solids4Foam ***********/

    label  nWedgePatches(0);
    vector wedgeDirVec(vector::zero);

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());
        }
    }

    reduce(nWedgePatches, maxOp<label>());

    if (nWedgePatches)
    {
        Info<< nl << "Axisymmetric case: disabling the solution in the "
            << "out-of-plane direction" << endl;

        // We will use const_cast to disable the out-of-plane direction
        Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                solD[cmpt] = -1;

                wordList dirs(3);
                dirs[0] = "x";
                dirs[1] = "y";
                dirs[2] = "z";
                Info<< "    out-of-plane direction: " << dirs[cmpt] << nl
                    << endl;
            }
        }
    }

    // Check all the face normals are in the same direction on the wedge patches
    // This is to avoid the case where a wedge patch is composed of two
    // disconnected regions with one on the front and one on the back
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            // Unit face normals on processor
            const vectorField nf = mesh.boundaryMesh()[patchI].faceNormals();

            if (nf.size() == 0)
            {
                FatalErrorIn("void Foam::mechanicsSubSolver::checkWedges() const")
                    << "There are no faces on the wedge patch "
                    << mesh.boundaryMesh()[patchI].name() << " on this processor:"
                    << nl << "Every processor should have at least one face on "
                    << "each wedge patch"
                    << abort(FatalError);
            }

            // Check that all the wedge face normals point in the same direction

            vector firstFaceNOnMasterProc = vector::zero;

            if (Pstream::master())
            {
                firstFaceNOnMasterProc = nf[0];
            }

            // Sync in parallel so that all processors have the master vector
            reduce(firstFaceNOnMasterProc, sumOp<vector>());

            forAll(nf, faceI)
            {
                if ((nf[faceI] & firstFaceNOnMasterProc) < 0)
                {
                    FatalErrorIn("void Foam::mechanicsSubSolver::checkWedges() const")
                        << "On wedge patch "
                        << mesh.boundaryMesh()[patchI].name() << " there are at "
                        << "least two faces with unit normals in the opposite "
                        << "directions" << nl
                        << "Please check that the wedge patches are correctly "
                        << "defined"
                        << abort(FatalError);
                }
            }
        }
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::mechanicsSubSolver::findRefCell()
{
    const vectorField& C = mesh_.C().internalField();
    const scalarField& V = mesh_.V().field();
    
    point C0 = gSum(C*V)/gSum(V);
    
    refCell_ = mesh_.findCell(C0);
    
    label count = (refCell_ >= 0);
    reduce(count, sumOp<label>());
    
    if (Pstream::master())
    {
        if (count == 0)
        {
            WarningIn("Foam::smallStrainSolver::findRefCell()")
                << "Failure to find reference cell." << endl;
        }
        else if (count > 1)
        {
            WarningIn("Foam::smallStrainSolver::findRefCell()")
                << "Reference cells found on multiple MPI tasks."
                << endl;
        }
    }
    
    if (count != 1)
    {
        refCell_ = -1;
    }
}


int Foam::mechanicsSubSolver::nCorrectors(const word fieldName)
{
    const dictionary& stressControl =
        mesh_.solutionDict().subDict("stressAnalysis");

    int nCorrectors = 1;
    if(stressControl.found("nCorrectors"))
    {
    #ifdef OPENFOAMFOUNDATION
        if(stressControl.isDict("nCorrectors"))
        {
            if (stressControl.subDict("nCorrectors").found(fieldName))
            {
                nCorrectors = stressControl.subDict("nCorrectors").lookup<int>(fieldName);
            }
            else if (stressControl.subDict("nCorrectors").found("default"))
            {
                nCorrectors = stressControl.subDict("nCorrectors").lookup<int>("default");
            }
            else
            {
                FatalErrorInFunction()
                << "Number of inner correctors missing for field \"" << fieldName << "\""
                << " in the fvSolution/stressAnalysis/nCorrectors sub-dict." << nl
                << "Please set this value using either the \"" << fieldName << "\"" 
                << " or \"default\" keywords" << nl << exit(FatalError);
            }
        }
        else
        {
            nCorrectors = stressControl.lookup<int>("nCorrectors");
        }
    #elif OPENFOAMESI      
        if(stressControl.isDict("nCorrectors"))
        {
            if (stressControl.subDict("nCorrectors").found(fieldName))
            {
                nCorrectors = stressControl.subDict("nCorrectors").get<int>(fieldName);
            }
            else if (stressControl.subDict("nCorrectors").found("default"))
            {
                nCorrectors = stressControl.subDict("nCorrectors").get<int>("default");
            }
            else
            {
                FatalErrorInFunction()
                << "Number of inner correctors missing for field \"" << fieldName << "\""
                << " in the fvSolution/stressAnalysis/nCorrectors sub-dict." << nl
                << "Please set this value using either the \"" << fieldName << "\"" 
                << " or \"default\" keywords" << nl << exit(FatalError);
            }
        }
        else
        {
            nCorrectors = stressControl.get<int>("nCorrectors");
        }
    #endif
    }

    return nCorrectors;
}


void Foam::mechanicsSubSolver::forceSummary()
{
    if (forceSummary_)
    {
        Info << "Force summary on domain boundaries:" << endl;

        scalar maxForce(0);
        vector totalForce(Foam::vector(0, 0, 0));

        forAll(mesh_.boundaryMesh(), patchID)
        {
            if(!isType<processorPolyPatch>(mesh_.boundaryMesh()[patchID]))
            {
                vector force;
                
                if(isType<tractionDisplacementFvPatchVectorField>(
                        D_.boundaryField()[patchID]))
                {
                    const tractionDisplacementFvPatchVectorField& pf
                        = refCast<const tractionDisplacementFvPatchVectorField>
                    (
                        D_.boundaryField()[patchID]
                    );
                    
                    force = pf.totalForce();
                }
                else
                {                
                    const symmTensorField& sigmaP = 
                        sigma_.boundaryField()[patchID];

                    const vectorField& patchSf = 
                        sigma_.boundaryField()[patchID].patch().Sf();

                    force = gSum(sigmaP & patchSf);
                }
            
                totalForce += force;
                maxForce = max(maxForce, mag(force));

                Info<< tab << mesh_.boundaryMesh()[patchID].name() 
                    << " (N) : " << force << endl;
            }
        }

        Info<< nl << tab << "Largest force (N) : " << maxForce << nl
            << tab << "Total force (N) : " << mag(totalForce) << nl
            << tab << "Error (-) : " << mag(totalForce)/maxForce << nl 
            << endl;
    }
}


void Foam::mechanicsSubSolver::checkWedges() const
{
    const fvMesh& mesh = this->mesh();

    label  nWedgePatches(0);
    vector wedgeDirVec(vector::zero);

    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            const wedgePolyPatch& wpp = refCast<const wedgePolyPatch>
            (
                mesh.boundaryMesh()[patchI]
            );

            nWedgePatches++;
            wedgeDirVec += cmptMag(wpp.centreNormal());
        }
    }

    reduce(nWedgePatches, maxOp<label>());

    if (nWedgePatches)
    {
        Info<< nl << "Axisymmetric case: disabling the solution in the "
            << "out-of-plane direction" << endl;

        // We will use const_cast to disable the out-of-plane direction
        Vector<label>& solD = const_cast<Vector<label>&>(mesh.solutionD());

        reduce(wedgeDirVec, sumOp<vector>());

        wedgeDirVec /= mag(wedgeDirVec);

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (wedgeDirVec[cmpt] > 1e-6)
            {
                solD[cmpt] = -1;

                wordList dirs(3);
                dirs[0] = "x";
                dirs[1] = "y";
                dirs[2] = "z";
                Info<< "    out-of-plane direction: " << dirs[cmpt] << nl
                    << endl;
            }
        }
    }

    // Check all the face normals are in the same direction on the wedge patches
    // This is to avoid the case where a wedge patch is composed of two
    // disconnected regions with one on the front and one on the back
    forAll(mesh.boundaryMesh(), patchI)
    {
        if (isA<wedgePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            // Unit face normals on processor
            const vectorField nf = mesh.boundaryMesh()[patchI].faceNormals();

            if (nf.size() == 0)
            {
                FatalErrorIn("void Foam::mechanicsSubSolver::checkWedges() const")
                    << "There are no faces on the wedge patch "
                    << mesh.boundaryMesh()[patchI].name() << " on this processor:"
                    << nl << "Every processor should have at least one face on "
                    << "each wedge patch"
                    << abort(FatalError);
            }

            // Check that all the wedge face normals point in the same direction

            vector firstFaceNOnMasterProc = vector::zero;

            if (Pstream::master())
            {
                firstFaceNOnMasterProc = nf[0];
            }

            // Sync in parallel so that all processors have the master vector
            reduce(firstFaceNOnMasterProc, sumOp<vector>());

            forAll(nf, faceI)
            {
                if ((nf[faceI] & firstFaceNOnMasterProc) < 0)
                {
                    FatalErrorIn("void Foam::mechanicsSubSolver::checkWedges() const")
                        << "On wedge patch "
                        << mesh.boundaryMesh()[patchI].name() << " there are at "
                        << "least two faces with unit normals in the opposite "
                        << "directions" << nl
                        << "Please check that the wedge patches are correctly "
                        << "defined"
                        << abort(FatalError);
                }
            }
        }
    }
}


Foam::tmp<Foam::volScalarField>
Foam::mechanicsSubSolver::calculateInitialGapWidth() const
{
    tmp<volScalarField> tgapWidth
    (
        new volScalarField
        (
            IOobject
            (
                "initialGapWidth",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar("initialGapWidth", dimless, GREAT)
        )
    );

    //- If required, update the AMI addressing
    forAll(D_.boundaryField(), patchI)
    {
        if
        (          
            isType<regionCoupledOFFBEATFvPatch>           
            (D_.boundaryField()[patchI].patch())
        )
        {        
            const regionCoupledOFFBEATFvPatch& patch
            = refCast<const regionCoupledOFFBEATFvPatch>         
            (D_.boundaryField()[patchI].patch());

            if(patch.owner())
            {
                patch.updateAMI();
            }
        }
    }

    //- Calculate the gap width on the slave
    forAll(D_.boundaryField(), patchI)
    {
        if
        (
            isType<regionCoupledOFFBEATFvPatch>
            (
                D_.boundaryField()[patchI].patch()
            )
        )
        {     
            const regionCoupledOFFBEATFvPatch& patch
                = refCast<const regionCoupledOFFBEATFvPatch>  
                (
                    D_.boundaryField()[patchI].patch()
                );

            if(!patch.owner())
            {         
                const regionCoupledOFFBEATFvPatch& nbrPatch
                    = refCast<const regionCoupledOFFBEATFvPatch>        
                (
                    patch.nbrPatch()
                );    
                    
                vectorField nf = patch.Sf() / patch.magSf();

                vectorField Cf = patch.Cf();
                               
                vectorField nbrCf =             
                patch.regionCoupledPatch().interpolate              
                (
                   nbrPatch.Cf()
                );

                scalarField& gWp = tgapWidth->boundaryFieldRef()[patchI];
                gWp = (nbrCf - Cf) & nf;
            }
        }
    }    

    //- Map gapWidth on the master side
    forAll(D_.boundaryField(), patchI)
    {
        if
        (   
            isType<regionCoupledOFFBEATFvPatch> 
            (
                D_.boundaryField()[patchI].patch()
            )
        )
        {    
            const regionCoupledOFFBEATFvPatch& patch
                = refCast<const regionCoupledOFFBEATFvPatch> 
                (
                    D_.boundaryField()[patchI].patch()
                );

            if(patch.owner())
            {           
                const regionCoupledOFFBEATFvPatch& nbrPatch
                    = refCast<const regionCoupledOFFBEATFvPatch> 
                    (
                        patch.nbrPatch()
                        );    

                const word nbrPatchName = 
                nbrPatch.name();
                const label nbrPatchID = 
                mesh_.boundaryMesh().findPatchID(nbrPatchName);

                scalarField& gWp = tgapWidth->boundaryFieldRef()[patchI];
                gWp =                    
                    patch.regionCoupledPatch().interpolate               
                    (
                        tgapWidth->boundaryField()[nbrPatchID]
                    );
            }               
        }
    }

    return tgapWidth;
}


void Foam::mechanicsSubSolver::correctCylindrical()
{
    if (cylindricalStress_)
    {
        // Create fields if not already available
        if (!sigmaCyl_.valid())
        {
            sigmaCyl_.set
            (
                new volSymmTensorField
                (
                    IOobject
                    (
                        "sigmaCyl",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedSymmTensor("sigmaCyl", dimPressure, symmTensor::zero)
                )
            );
        }

        if (!epsilonCyl_.valid())
        {
            epsilonCyl_.set
            (
                new volSymmTensorField
                (
                    IOobject
                    (
                        "epsilonCyl",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedSymmTensor("epsilonCyl", dimPressure, symmTensor::zero)
                )
            );
        }
        
        if (!DCyl_.valid())
        {
            DCyl_.set
            (
                new volVectorField
                (
                    IOobject
                    (
                        "DCyl",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("DCyl", dimLength, vector::zero)
                )
            );
        }
        
        // Working fields
        vectorField& DCyl = DCyl_->ref();
        symmTensorField& sigmaCyl = sigmaCyl_->ref();
        symmTensorField& epsilonCyl = epsilonCyl_->ref();
        
        const vectorField& Di = D_.internalField();
        const symmTensorField& sigmai = sigma_.internalField();
        const symmTensorField& epsiloni = epsilon_.internalField();
        const vectorField& C = mesh_.C().internalField();
        
        // Hardcoded for rotation around the z-axis
        // TODO - Needs to be generalised for arbitrary rotation axis
        forAll(sigmai, addrI)
        {
            tensor R = tensor::zero;

            const vector& Ci = C[addrI];
            const scalar radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1]);

            R.xx() = Ci[0]/radius;
            R.xy() = Ci[1]/radius;
            R.yy() = Ci[0]/radius;
            R.yx() = -Ci[1]/radius;
            R.zz() = 1;

            DCyl[addrI] = (R & Di[addrI]);
            sigmaCyl[addrI] = symm(R & sigmai[addrI] & R.T());
            epsilonCyl[addrI] = symm(R & epsiloni[addrI] & R.T());
        }

        forAll(sigma_.boundaryField(), patchI)
        {
            vectorField& DCylP = DCyl_->boundaryFieldRef()[patchI];
            symmTensorField& sigmaCylP = sigmaCyl_->boundaryFieldRef()[patchI];
            symmTensorField& epsilonCylP = epsilonCyl_->boundaryFieldRef()[patchI];
            
            const vectorField& DP = D_.boundaryField()[patchI];
            const symmTensorField& sigmaP = sigma_.boundaryField()[patchI];
            const symmTensorField& epsilonP = epsilon_.boundaryField()[patchI];
            const fvPatch& p = sigma_.boundaryField()[patchI].patch();
            const vectorField& Cf = p.Cf();

            forAll(sigmaP, faceI)
            {
                tensor R = tensor::zero;

                const vector& Ci = Cf[faceI];
                const scalar radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1]);

                R.xx() = Ci[0]/radius;
                R.xy() = Ci[1]/radius;
                R.yy() = Ci[0]/radius;
                R.yx() = -Ci[1]/radius;
                R.zz() = 1;

                DCylP[faceI] = (R & DP[faceI]);
                sigmaCylP[faceI] = symm(R & sigmaP[faceI] & R.T());
                epsilonCylP[faceI] = symm(R & epsilonP[faceI] & R.T());
            }
        }

        DCyl_->correctBoundaryConditions();
        sigmaCyl_->correctBoundaryConditions();
        epsilonCyl_->correctBoundaryConditions();
    }
}

void Foam::mechanicsSubSolver::correctSpherical()
{
    if (sphericalStress_)
    {
        // Create fields if not already available
        if (!sigmaSph_.valid())
        {
            sigmaSph_.set
            (
                new volSymmTensorField
                (
                    IOobject
                    (
                        "sigmaSph",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedSymmTensor("sigmaSph", dimPressure, symmTensor::zero)
                )
            );
        }

        if (!epsilonSph_.valid())
        {
            epsilonSph_.set
            (
                new volSymmTensorField
                (
                    IOobject
                    (
                        "epsilonSph",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedSymmTensor("epsilonSph", dimPressure, symmTensor::zero)
                )
            );
        }

        if (!DSph_.valid())
        {
            DSph_.set
            (
                new volVectorField
                (
                    IOobject
                    (
                        "DSph",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    mesh_,
                    dimensionedVector("DSph", dimLength, vector::zero)
                )
            );
        }

        // Working fields
        vectorField& DSph = DSph_->ref();
        symmTensorField& sigmaSph = sigmaSph_->ref();
        symmTensorField& epsilonSph = epsilonSph_->ref();

        const vectorField& Di = D_.internalField();
        const symmTensorField& sigmai = sigma_.internalField();
        const symmTensorField& epsiloni = epsilon_.internalField();
        const vectorField& C = mesh_.C().internalField();

        // Hardcoded for rotation around the z-axis
        // TODO - Needs to be generalised for arbitrary rotation axis
        forAll(sigmai, addrI)
        {
            const vector& Ci = C[addrI];
            const scalar radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1] + Ci[2]*Ci[2]);

            scalar theta, phi;

            if (radius <= 1e-10) //- For the element where its center is the center of the sphere.
            {
              theta = 0.0;
              phi = 0.0;
            }
            else
            {
              theta = acos(Ci[2]/radius);
              phi = atan2(Ci[1],Ci[0]);
            }

            const tensor R(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta),
                           cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta),
                           -sin(phi),           cos(phi),            0.0);

            DSph[addrI] = (R & Di[addrI]);
            sigmaSph[addrI] = symm(R & sigmai[addrI] & R.T());
            epsilonSph[addrI] = symm(R & epsiloni[addrI] & R.T());
        }

        forAll(sigma_.boundaryField(), patchI)
        {
            vectorField& DSphP = DSph_->boundaryFieldRef()[patchI];
            symmTensorField& sigmaSphP = sigmaSph_->boundaryFieldRef()[patchI];
            symmTensorField& epsilonSphP = epsilonSph_->boundaryFieldRef()[patchI];

            const vectorField& DP = D_.boundaryField()[patchI];
            const symmTensorField& sigmaP = sigma_.boundaryField()[patchI];
            const symmTensorField& epsilonP = epsilon_.boundaryField()[patchI];
            const fvPatch& p = sigma_.boundaryField()[patchI].patch();
            const vectorField& Cf = p.Cf();

            forAll(sigmaP, faceI)
            {

                const vector& Ci = Cf[faceI];
                const scalar radius = sqrt(Ci[0]*Ci[0] + Ci[1]*Ci[1] + Ci[2]*Ci[2]);

                scalar theta, phi;

                if (radius <= SMALL) //- For the element where its center is the center of the sphere.
                {
                  theta = 0.0;
                  phi = 0.0;
                }
                else
                {
                  theta = acos(Ci[2]/radius);
                  phi = atan2(Ci[1],Ci[0]);
                }

                const tensor R(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta),
                               cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta),
                               -sin(phi),           cos(phi),            0.0);

                DSphP[faceI] = (R & DP[faceI]);
                sigmaSphP[faceI] = symm(R & sigmaP[faceI] & R.T());
                epsilonSphP[faceI] = symm(R & epsilonP[faceI] & R.T());
            }
        }

        DSph_->correctBoundaryConditions();
        sigmaSph_->correctBoundaryConditions();
        epsilonSph_->correctBoundaryConditions();
    }
}

Foam::tmp<Foam::surfaceScalarField>
Foam::mechanicsSubSolver::correctedTwoMuLambdaF
(
    const volScalarField& twoMuLambda
) const
{
    if (multiMaterialCorrection_)
    {
        return harmonic(mesh()).interpolate(twoMuLambda);
    }
    else
    {
        Istream& schemeData(mesh().laplacianScheme("laplacian(DD,D)"));

        // Advance Istream past laplacian name (Gauss) and move to 
        // interpolation scheme name (linear, harmonic etc.)
        const word schemeName(schemeData);

        tmp<surfaceInterpolationScheme<scalar>> 
        tinterpGammaScheme_ = tmp<surfaceInterpolationScheme<scalar>>
        (
            surfaceInterpolationScheme<scalar>::New(mesh(), schemeData)
        );

        return tinterpGammaScheme_().interpolate(twoMuLambda);
    }
}


Foam::tmp<Foam::volTensorField>
Foam::mechanicsSubSolver::correctedGradD
(
    const volScalarField& twoMuLambda,
    const volVectorField& D,
    const volTensorField& gradD,
    const volTensorField& sigmaExp
) const
{
    if (multiMaterialCorrection_)
    {
        // Multi-material interface weights (both for the gradient and strain
        // parts of sigmaExp)
        const surfaceScalarField& w_eps = materialInterface_->weightsEps();
        const surfaceScalarField& w_grad = materialInterface_->weightsGrad();

        // Mesh weights and normals
        const surfaceScalarField& beta = mesh().weights();
        surfaceVectorField nf = mesh().Sf()/mesh().magSf();

        if( mesh_.foundObject<volTensorField>("F") )
        {
            const volTensorField& Finv = 
            mesh_.foundObject<fvMesh>("referenceMesh") ?
            mesh_.lookupObject<volTensorField>("relFinv") :
            mesh_.lookupObject<volTensorField>("Finv");

            const volScalarField& J = 
            mesh_.foundObject<fvMesh>("referenceMesh") ?
            mesh_.lookupObject<volScalarField>("relJ") :
            mesh_.lookupObject<volScalarField>("J");

            surfaceTensorField Finvf(fvc::interpolate(Finv));
            surfaceScalarField Jf(fvc::interpolate(J));

            nf = Jf*Finvf.T() & nf;
            nf /= mag(nf);
        }

        surfaceInterpolation interp(mesh());
        const surfaceScalarField& rDelta = interp.nonOrthDeltaCoeffs();

        // Delta vectors
        surfaceVectorField deltaf = mesh().delta();

        // Non-orthogonal correction vectors
        surfaceVectorField k = (nf - deltaf/(nf&deltaf));
        
        // Split sigmaExp into gradient and (in-elastic) strain parts
        const volScalarField& mu(mat_.mu());
        const volScalarField& lambda(mat_.lambda());
        volTensorField sigmaExpGrad = mu*gradD.T() - (mu + lambda)*gradD + lambda*tr(gradD)*I;
        volTensorField sigmaExpEps = sigmaExp - sigmaExpGrad;

        // Interpolate the displacement field to the faces and apply
        // multi-material interface corrections
        surfaceVectorField Df
        (
            D.name(), 
            ( 
                fvc::interpolate(twoMuLambda*D, "interpolate(" + D.name() + ")")
                + 
                w_eps*beta*(1-beta)/(rDelta*rDelta)*
                (
                    (nf & fv::uncorrectedSnGrad<tensor>(mesh()).snGrad(sigmaExpEps))
                )
                +
                w_grad*beta*(1-beta)/(rDelta*rDelta)*
                (
                    (nf & fv::uncorrectedSnGrad<tensor>(mesh()).snGrad(sigmaExpGrad)) 
                    + (k & fv::uncorrectedSnGrad<tensor>(mesh()).snGrad
                        ((twoMuLambda*gradD)()))
                )
            ) / fvc::interpolate(twoMuLambda, "interpolate(" + D.name() + ")") 
        );

        // Calculate gradient of displacement
        
        tmp<volTensorField> tGrad
        (
            new volTensorField
            (
                IOobject
                (
                    "tGrad",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                // fvc::grad(Df) 
                fvc::leastSquaresGradf(Df, D, "grad(D" + D.name() + ")")
            )
        );

        volTensorField& tGradRef(tGrad.ref());

        // fv::gaussGrad<vector>::correctBoundaryConditions(D, tGradRef);

        // Force x-y gradient for 1.5D simulations
        // The slice-wise constant epsilonZZ value  is then added by the 
        // correctModifiedPlaneStrain() function in rheology
        if(rheo_.modifiedPlaneStrain())
        {
            tensor planeStrainTensor = Foam::tensor(1, 0, 0, 0, 1, 0, 0, 0, 0);
            tGradRef = (planeStrainTensor & (tGradRef & planeStrainTensor));
        }

        return tGrad;
    }
    else
    {
        return fvc::grad(D, "grad(D" + D.name() + ")");
    }
}


Foam::tmp<Foam::volTensorField>
Foam::mechanicsSubSolver::correctedGradDD
(
    const volScalarField& twoMuLambda,
    const volVectorField& DD
) const
{
    if (multiMaterialCorrection_)
    {
        // Interpolate the displacement field to the faces and apply
        // multi-material interface corrections
        surfaceVectorField DDf("DD", (
            fvc::interpolate(twoMuLambda*DD, "interpolate(DD)")
        ) / fvc::interpolate(twoMuLambda, "interpolate(DD)"));
        
        tmp<volTensorField> tGrad
        (
            new volTensorField
            (
                IOobject
                (
                    "tGrad",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                fvc::grad(DDf)
            )
        );

        volTensorField& tGradRef(tGrad.ref());

        fv::gaussGrad<vector>::correctBoundaryConditions(DD, tGradRef);

        // Calculate gradient of displacement
        return tGrad;
    }
    else
    {
        return fvc::grad(DD, "grad(DD)");
    }
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::mechanicsSubSolver::fvmLaplacianD
(
    const surfaceScalarField& twoMuLambda_f,
    const volVectorField& D,
    const word& laplacianName
) const
{
    if(rheo_.modifiedPlaneStrain())
    {
        surfaceTensorField twoMuLambda_f_tensor = 
        twoMuLambda_f*Foam::tensor(1, 0, 0, 0, 1, 0, 0, 0, 0);

        return fvmLaplacianD(twoMuLambda_f_tensor, D, laplacianName);
    }

    // In case of multi-material, impose uncorrected Gauss laplacian scheme.
    // Note that the surface-interpolation scheme is not used as twoMuLambda
    // is already provided as a surfaceScalarField.
    if (multiMaterialCorrection_)
    {
        tmp<surfaceInterpolationScheme<scalar>> 
        tinterpGammaScheme(new harmonic(mesh()));

        tmp<fv::snGradScheme<vector>> 
        tsnGradScheme(new fv::uncorrectedSnGrad<vector>(mesh()));

        return fv::gaussLaplacianScheme<vector, scalar>(
            mesh(), tinterpGammaScheme, tsnGradScheme).fvmLaplacian
            (twoMuLambda_f, D);
    }
    else
    {
        return fvm::laplacian(twoMuLambda_f, D, laplacianName);
    }
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::mechanicsSubSolver::fvmLaplacianD
(
    const surfaceTensorField& twoMuLambda_f,
    const volVectorField& D,
    const word& laplacianName
) const
{
    // In case of multi-material, impose uncorrected Gauss laplacian scheme.
    // Note that the surface-interpolation scheme is not used as twoMuLambda
    // is already provided as a surfaceScalarField.
    if (multiMaterialCorrection_)
    {
        tmp<surfaceInterpolationScheme<tensor>> 
        tinterpGammaScheme(new linear<tensor>(mesh()));

        tmp<fv::snGradScheme<vector>> 
        tsnGradScheme(new fv::uncorrectedSnGrad<vector>(mesh()));

        return fv::gaussLaplacianScheme<vector, tensor>(
            mesh(), tinterpGammaScheme, tsnGradScheme).fvmLaplacian
            (twoMuLambda_f, D);
    }
    else
    {
        Istream& schemeData(mesh().laplacianScheme("laplacian(DD,D)"));

        // Advance Istream past laplacian name (Gauss) and 
        // interpolation scheme name
        const word schemeName(schemeData);
        const word interpSchemeName(schemeData);
            
        // As the surface-interpolation scheme is not used (twoMuLambda
        // is already provided as a surfaceScalarField) we force it to be
        // linear to avoid problems when the users pick an interpolation
        // scheme that is not compatible with a tensorial laplacian coeff.
        tmp<surfaceInterpolationScheme<tensor>> 
        tinterpGammaScheme(new linear<tensor>(mesh()));

        tmp<fv::snGradScheme<vector>> 
        tsnGradScheme(fv::snGradScheme<vector>::New(mesh(), schemeData));
        
        return fv::gaussLaplacianScheme<vector, tensor>(
            mesh(), tinterpGammaScheme, tsnGradScheme).fvmLaplacian
            (twoMuLambda_f, D);
    }
}


template<>
Foam::tmp<Foam::volTensorField>
Foam::mechanicsSubSolver::correctedSigmaExp
(
    const volTensorField& sigma,
    const volScalarField& twoMuLambda,
    const volTensorField& gradD
) const;


template<>
Foam::tmp<Foam::volTensorField>
Foam::mechanicsSubSolver::correctedSigmaExp
(
    const volSymmTensorField& sigma,
    const volScalarField& twoMuLambda,
    const volTensorField& gradD
) const;


template<>
Foam::tmp<Foam::volVectorField>
Foam::mechanicsSubSolver::correctedDivSigmaExp
(
    const volTensorField& sigmaExp,
    const volSymmTensorField& sigma,
    const volScalarField& twoMuLambda,
    const surfaceScalarField& twoMuLambda_f,
    const volVectorField& D,
    const volTensorField& gradD
) const;


template<>
Foam::tmp<Foam::volVectorField>
Foam::mechanicsSubSolver::correctedDivSigmaExp
(
    const volTensorField& sigmaExp,
    const volTensorField& sigma,
    const volScalarField& twoMuLambda,
    const surfaceScalarField& twoMuLambda_f,
    const volVectorField& D,
    const volTensorField& gradD
) const;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mechanicsSubSolver::mechanicsSubSolver
(
    fvMesh& mesh,
    const materials& materials,
    rheology& rheo,
    const dictionary& mechanicsDict
)
:     
    mesh_(mesh),
    mat_(materials),
    rheo_(rheo),
    mechanicsDict_(mechanicsDict),
    D_
    (
        IOobject
        (
            "D",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("D", dimLength, vector::zero)
    ),
    gradD_
    (
        IOobject
        (
            "gradD",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("gradD", dimless, tensor::zero)
    ),
    sigma_
    (
        IOobject
        (
            "sigma",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("sigma", dimPressure, symmTensor::zero)
    ),
    sigmaExp_
    (
        IOobject
        (
            "sigmaExp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedTensor("sigmaExp", dimPressure, tensor::zero)
    ),
    divSigmaExp_
    (
        IOobject
        (
            "divSigmaExp",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedVector("divSigmaExp", dimPressure/dimLength, vector::zero)
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("epsilon", dimless, symmTensor::zero)
    ),
    interfaceP_
    (
        IOobject
        (
            "interfaceP",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("interfaceP", dimPressure, 0)
    ),
    gapWidth_
    (
        IOobject
        (
            "gapWidth",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        calculateInitialGapWidth()
    ),
    RhieChowCorrection_(true),
    residual_(vector::zero),
    initialResidual_(vector::zero),
    relResidual_(0),
    initialRelResidual_(0),
    absErr_(0.0),
    initialAbsErr_(0.0),
    refCell_(-1),
    multiMaterialCorrection_(false),
    forceSummary_(false),
    cylindricalStress_(false),
    sphericalStress_(false)
{
    // If the case is axisymmetric, disable solving in the out-of-plane
    checkWedges();
    
    if
    (
#ifdef OPENFOAMFOUNDATION    
        !mechanicsDict_.isNull()
#elif OPENFOAMESI
        !mechanicsDict_.isNullDict()
#endif
    )
    {
        forceSummary_ = 
            mechanicsDict_.lookupOrDefault<bool>("forceSummary", false);
        
        cylindricalStress_ = 
            mechanicsDict_.lookupOrDefault<bool>("cylindricalStress", false);

        sphericalStress_ =
            mechanicsDict_.lookupOrDefault<bool>("sphericalStress", false);

        RhieChowCorrection_ =
            mechanicsDict_.lookupOrDefault<bool>("RhieChowCorrection", true);
        
        if (mechanicsDict_.found("multiMaterialCorrection"))
        {
            multiMaterialCorrection_ = true;
            
            materialInterface_ = multiMaterialInterface::New
            (
                mesh,
                mechanicsDict_.subDict("multiMaterialCorrection")
            );
        }
    }
}


Foam::autoPtr<Foam::mechanicsSubSolver>
Foam::mechanicsSubSolver::New
(
    fvMesh& mesh,
    const materials& materials,
    rheology& rheology,
    const dictionary& solverDict
)
{
    // Default pick type
    word type("constant");

    // Prepare options dict
    dictionary mechanicsDict
    (
        solverDict.subOrEmptyDict("mechanicsSolverOptions")
    );

    if (solverDict.found("mechanicsSolver"))
    {
        if (solverDict.isDict("mechanicsSolver"))
        {
            // Case 1: mechanicsSolver is a dictionary
            mechanicsDict = solverDict.subDict("mechanicsSolver");
            mechanicsDict.lookup("type") >> type;
        }
        else
        {
            // Case 2: mechanicsSolver is a word (type)
            solverDict.lookup("mechanicsSolver") >> type;
            // mechanicsDict stays as mechanicsSolverOptions
        }
    }

    // Map deprecated `fromLatestTime` into new name `constant`
    if (type == "fromLatestTime")
    {
        WarningIn("mechanicsSubSolver::New(const fvMesh&, const materials&, rheology&, const dictinoary&)")
            << "Type 'fromLatestTime' for `mechanicsSolver` is deprecated. " 
            << "Please use 'constant' instead." << nl << endl;
        
        type = "constant";
    }
    
    Info << "Selecting mechanicsSolver model " << type << endl;
    
    auto cstrIter
        = dictionaryConstructorTablePtr_->find(type);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "mechanicsSubSolver::New(const fvMesh&, const materials&," 
            "const rheology&, const dictionary&)"
        )
            << "Unknown mechanicsSubSolver type "
            << type << nl << nl
            << "Valid types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<mechanicsSubSolver>
    (
        cstrIter()
        (
            mesh, materials, rheology, mechanicsDict
        )
    );
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mechanicsSubSolver::~mechanicsSubSolver()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


void Foam::mechanicsSubSolver::updateTotalFields()
{    
    rheo_.updateTotalFields();
}


void Foam::mechanicsSubSolver::writeStressFields()
{
    //- Calculate the von Mises stress
    volScalarField sigmaEq
    (
        IOobject
        (
            "sigmaEq",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        sqrt((3.0/2.0)*magSqr(dev(sigma_)))
    );
    
    sigmaEq.write();

    Info<< "Max von Mises stress (sigmaEq) = " << max(sigmaEq).value()
        << endl;
}


bool Foam::mechanicsSubSolver::converged
(
    volVectorField& D, 
    const vector& residual,
    const scalar& relResidual,
    const scalar& absErr
)
{
    bool converged(false);

    const dictionary& stressControl = 
    mesh_.solutionDict().subDict("stressAnalysis");

    bool useRelRes(stressControl.lookupOrDefault("useRelResD", true));
    
    vector convergenceTolerance
    (
        stressControl.lookupOrDefault(D.name(), vector(1e-6,1e-6,1e-6))
    );
    
    scalar absErrTolerance
    (
        stressControl.lookupOrDefault("absErrD", 0.0)
    );

    if(absErr > absErrTolerance)
    {
        if ( useRelRes )
        {
            scalar relConvergenceTolerance
            (
                readScalar(stressControl.lookup("relD"))
            );

            // Consider converged if both initialResidual and relResidual are 
            // below threshold or if one of the two is 10 times below threshold.
            converged = 
            ( 
                residual < convergenceTolerance 
                && 
                relResidual < relConvergenceTolerance 
            )
            or
            (
                residual < convergenceTolerance/10
            )
            or
            (
                relResidual < relConvergenceTolerance/10 
            )
            ; 
        }
        else
        {
            // In this case only the residual computed by OpenFOAM is accounted for 
            // the convergence monitoring.
            converged = 
            ( 
                residual < convergenceTolerance 
            ); 
        }
    }
    else
    {
        converged = true;
    }

    return converged;
}


bool Foam::mechanicsSubSolver::converged(volVectorField& D)
{
    return converged(D, initialResidual_, initialRelResidual_, initialAbsErr_);
}


void Foam::mechanicsSubSolver::coupleBoundaries
(
    volVectorField& D, 
    bool coupleOrUnCouple
)
{
    forAll(D.boundaryField(), patchI)
    {
        if (isA<implicitGapContactFvPatchVectorField>(
            D.boundaryField()[patchI]))
        {
            refCast<implicitGapContactFvPatchVectorField>(
                D.boundaryFieldRef()[patchI]).coupleField(coupleOrUnCouple);
        }
    }
}

// ************************************************************************* //

