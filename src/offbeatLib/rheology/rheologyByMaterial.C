// ************************************************************************* //


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

#include "rheologyByMaterial.H"
#include "zeroGradientFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "emptyFvPatchFields.H"
#include "gapGasModel.H"

#include "Time.H"
#include "regionCoupledOFFBEATFvPatch.H"
#include "PrimitivePatchInterpolation.H"

#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rheologyByMaterial, 0);
    addToRunTimeSelectionTable
    (
        rheology, 
        rheologyByMaterial, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rheologyByMaterial::rheologyByMaterial
(
    const fvMesh& mesh,
    materials& mat,
    const dictionary& rheologyDict
)
:
    rheology(mesh, mat, rheologyDict),
    precision_(0),
    springModulus_(0),
    springPreCompression_(0),
    mapper_(nullptr),
    frictionCoeff_(0),
    penaltyScaleFactFriction_(0.001),
    frictionSaturationPressure_(1e6),
    timeCheck_(0),
    relax_(1.0),
    springForce_(),
    pressureForce_(),
    frictionForce_(),
    slip0_(),
    slip_()
{    
    lawsList_.clear();
    lawsList_.resize(mat_.materialsList().size());

    if(modifiedPlaneStrain_)
    {
        precision_ = 
        rheologyDict_.lookupOrDefault<scalar>("precisionSpring", 1e-6);

        springModulus_ = 
        rheologyDict_.lookupOrDefault<scalar>("springModulus", 3.5e3);  

        springPreCompression_ = 
        rheologyDict_.lookupOrDefault<scalar>("springPreCompression", 0.0);  

        frictionCoeff_ = 
        rheologyDict_.lookupOrDefault<scalar>("frictionCoefficient", 0.0); 

        penaltyScaleFactFriction_ = 
        rheologyDict_.lookupOrDefault<scalar>("penaltyFriction", 0.001);

        frictionSaturationPressure_ = 
        rheologyDict_.lookupOrDefault<scalar>("frictionSaturationPressure", 1e6);

        relax_ = 
        rheologyDict_.lookupOrDefault<scalar>("relax", 1.0);


#ifdef OPENFOAMFOUNDATION        
        // Read fluid pressure list
        fluidPressureList_.set
        ( 
            new Function1s::Table<scalar>
            (
                "coolantPressureList", rheologyDict_.subDict("coolantPressureList")
            )
        );    
#elif OPENFOAMESI
        // Read fluid pressure list
        fluidPressureList_.reset
        ( 
            Function1<scalar>::New
            (
                "table", rheologyDict_.subDict("coolantPressureList"), "table"
            )
        ); 
#endif            
    }
    
    forAll(mat_.materialsList(), i)
    {
        const materialModel& materialI(mat_.materialsList()[i]);

        const word& matName(materialI.name());
        const dictionary& materialModelDict(materialI.materialModelDict());
        
        Info << "     For material " << matName << ": " << endl;
        
        lawsList_.set
        (
            i,
            constitutiveLaw::New
            (
                mesh_,
                materialModelDict
            )
        );

        Info << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::rheologyByMaterial::~rheologyByMaterial()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::rheologyByMaterial::correct
(
    volSymmTensorField& sigma, 
    volSymmTensorField& epsilon,
    volVectorField& D
)
{        
    const labelListList& matAddrList(mat_.matAddrList());

    if(thermalExpansion_)
    {
        epsTh_ = mat_.thermalExpansion();

        if(mesh_.foundObject<surfaceScalarField>("multiMaterialWeightsEps"))
        {
            surfaceScalarField w = mesh_.lookupObject<surfaceScalarField>("multiMaterialWeightsEps");

            forAll(epsTh_.boundaryField(), patchi)
            {
                if (!epsTh_.boundaryField()[patchi].coupled())
                {
                    const labelList faceCells = mesh_.boundary()[patchi].faceCells();
                    symmTensorField epsThi(epsTh_.internalField(), faceCells);
                    const scalarField wb = w.boundaryField()[patchi];
                    symmTensorField& epsThb = epsTh_.boundaryFieldRef()[patchi];
                    epsThb = (1-wb)*epsThb + wb*epsThi;
                }
            }    
        }
    }


    // For planeStress or modifiedPlaneStrain, correct strain in the out of 
    // plane direction
    if (planeStress_)
    {
        if (mesh_.solutionD()[vector::Z] > -1)
        {            
            FatalErrorInFunction()
                << "For planeStress, this material law assumes the empty "
                << "direction is the Z direction!" << exit(FatalError);
        }

        correctPlaneStress(epsilon);
    }
    else if(modifiedPlaneStrain_)
    {        
        if (mesh_.solutionD()[vector::Z] > -1)
        {
            
            FatalErrorInFunction()
                << "For modifiedPlaneStrain, this material law assumes the empty "
                << "direction is the Z direction!" << exit(FatalError);
        }

        correctModifiedPlaneStrain(epsilon, D);
    }


    if(mesh_.foundObject<surfaceScalarField>("multiMaterialWeightsEps"))
    {
        surfaceScalarField w = mesh_.lookupObject<surfaceScalarField>("multiMaterialWeightsEps");

        forAll(additionalStrain_.boundaryField(), patchi)
        {
            if (!additionalStrain_.boundaryField()[patchi].coupled())
            {
                const labelList faceCells = mesh_.boundary()[patchi].faceCells();
                symmTensorField additionalStrain_i(additionalStrain_.internalField(), faceCells);
                const scalarField wb = w.boundaryField()[patchi];
                symmTensorField& additionalStrain_b = additionalStrain_.boundaryFieldRef()[patchi];
                additionalStrain_b = (1-wb)*additionalStrain_b + wb*additionalStrain_i;
            }
        }    
    }
    // Set elastic strain tensor
    epsEl_ = epsilon - additionalStrain_;

    // Subtract thermal strain
    if(thermalExpansion_)
    {
        epsEl_ -= epsTh_;
    }  

    volScalarField sigmaHydExplicit(sigmaHyd_);
    
    // Correct the elastic strain tesnor removing oldTime copmonents of plastic
    // strains
    forAll(mat_.materialsList(), i)
    {       
        const labelList& addr(matAddrList[i]);
        lawsList_[i].correctEpsilonEl(epsEl_, addr);
    }
    
    // Correct material rheology laws. It calculates the hydrostatic and deviatoric
    // component of the stress and it updates the elastic strain component.
    // Update damage models where nedded.
    forAll(mat_.materialsList(), i)
    {        
        const labelList& addr(matAddrList[i]);
        lawsList_[i].correct(sigmaHydExplicit, sigmaDev_, epsEl_, addr);
        mat_.materialsListRef()[i].correctDamage(sigmaHydExplicit, sigmaDev_, epsEl_);
    }

    if(solvePressureEqn_)
    {
        solvePressureEqn(sigmaHydExplicit);
    }
    else
    {
        sigmaHyd_ = sigmaHydExplicit;
    }

    // Calculate sigma from combination of sigmaHyd and sigmaDev.
    forAll(mat_.materialsList(), i)
    {
        const labelList& addr(matAddrList[i]);
        lawsList_[i].updateStress(sigma, sigmaHyd_, sigmaDev_, addr);
    }

    sigma.correctBoundaryConditions();
}


void Foam::rheologyByMaterial::updateTotalFields()
{

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        const labelList& addr(matAddrList[i]);
        lawsList_[i].updateTotalFields(addr);
    }       
}


Foam::volSymmTensorField Foam::rheologyByMaterial::getLawsAdditionalStrain()
{
    volSymmTensorField lawsAdditionalStrain
    (
        IOobject
        (
            "lawsAdditionalStrain",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedSymmTensor("lawsAdditionalStrain", dimless, symmTensor::zero),
        calculatedFvPatchField<scalar>::typeName
    );

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        const labelList& addr(matAddrList[i]);
        lawsList_[i].correctAdditionalStrain(lawsAdditionalStrain, addr);
    }

    return lawsAdditionalStrain;
}


void Foam::rheologyByMaterial::correctPlaneStress
(
    volSymmTensorField& epsilon
)
{
    volScalarField nu("nu", mat_.nu());

    // Get additional strains defined by the constitutive mechanical laws
    volSymmTensorField lawsAdditionalStrain
    (
        getLawsAdditionalStrain()
    );

    epsilon.replace
    (
        symmTensor::ZZ,
       -(nu/(1-nu))
       *(epsilon.component(symmTensor::XX) + epsilon.component(symmTensor::YY))
       +(nu/(1-nu))
       *
       (
            additionalStrain_.component(symmTensor::XX) +
            additionalStrain_.component(symmTensor::YY) + 
            lawsAdditionalStrain.component(symmTensor::XX) + 
            lawsAdditionalStrain.component(symmTensor::YY)
        )
       +(
            additionalStrain_.component(symmTensor::ZZ) +
            lawsAdditionalStrain.component(symmTensor::ZZ)
        )
    );

    if(thermalExpansion_)
    {
        epsilon.replace
        (
            symmTensor::ZZ,
            epsilon.component(symmTensor::ZZ)
           +(nu/(1-nu))
           *(epsTh_.component(symmTensor::XX) + epsTh_.component(symmTensor::YY))
           +(epsTh_.component(symmTensor::ZZ))
        );
    }

    // const labelListList& matAddrList(mat_.matAddrList());

    // forAll(mat_.materialsList(), i)
    // {
    //     const labelList& addr(matAddrList[i]);
    //     lawsList_[i].correctAxialStrain(epsilon, addr);
    // }
}


void Foam::rheologyByMaterial::correctModifiedPlaneStrain
(
    volSymmTensorField& epsilon,
    volVectorField& D
)
{
    // Modified plane strain (slice-based axial equilibrium)
    //
    // We enforce axial force balance per slice (cellZone) by solving for a
    // slice-uniform axial strain epsZ (= epsilon.zz()) such that:
    //
    //   ∫_slice sigmaZ dA = F_spring + F_pressure + F_friction
    //
    // With isotropic linear elasticity and “additional strains” (thermal, inelastic,
    // etc.), the axial stress can be written in the form:
    //
    //   sigmaZ = (2*mu + lambda)*epsZ  - (2*mu + lambda)*epsAddZ  + lambda*(epsElX + epsElY)
    //
    // Assuming epsZ is constant within each slice, we replace surface integrals by
    // slice averages times area. The slice area is obtained from geometry as:
    //
    //   A_slice = V_slice / (angularFraction * H_slice)
    //
    // This yields an explicit expression for epsZ_slice as a function of slice-averaged
    // material fields and external axial forces (spring / gas / friction). Because
    // forces depend on the axial displacement (and displacement depends on epsZ),
    // we iterate with relaxation until the average axial force converges.

    const scalar pi = Foam::constant::mathematical::pi;

    if (!mesh_.foundObject<sliceMapper>("sliceMapper"))
    {
        return;
    }

    if (mapper_ == nullptr)
    {
        mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
    }

    // Gap gas pressure contribution (top/bottom force balance)
    const gapGasModel& gapGas =
        mesh_.lookupObject<gapGasModel>("gapGas");

    const globalOptions& globalOpt =
        mesh_.lookupObject<globalOptions>("globalOptions");

    // Angular fraction used to scale slice quantities to the 360° rod
    const scalar angularFraction = globalOpt.angularFraction();
    
    // Pin direction
    const vector pinDirection = globalOpt.pinDirection();

    // Number of slices (equivalent to cellZones)
    const int nSlices = mapper_->nSlices();

    if (nSlices <= 0)
    {
        FatalErrorInFunction()
            << "When the option 'modifiedPlaneStrain' is activated, the number "
            << "of slices cannot be less than 1," << nl
            << "as the slices are needed to perform the axial force balance." << nl
            << "Check that a sliceMapper different than \"none\" is used "
            << "or check the slice definition."
            << abort(FatalError);
    }

    if (!springForce_.valid())
    {
#ifdef OPENFOAMFOUNDATION
        springForce_.set(new scalarField(nSlices, 0.0));
#elif OPENFOAMESI
        springForce_.reset(new scalarField(nSlices, 0.0));
#endif
    }

    if (!pressureForce_.valid())
    {
#ifdef OPENFOAMFOUNDATION
        pressureForce_.set(new scalarField(nSlices, 0.0));
#elif OPENFOAMESI
        pressureForce_.reset(new scalarField(nSlices, 0.0));
#endif
    }

    if (!frictionForce_.valid())
    {
#ifdef OPENFOAMFOUNDATION
        frictionForce_.set(new scalarField(nSlices, 0.0));
#elif OPENFOAMESI
        frictionForce_.reset(new scalarField(nSlices, 0.0));
#endif
    }

    if (!slip0_.valid())
    {
#ifdef OPENFOAMFOUNDATION
        slip0_.set(new scalarField(nSlices, 0.0));
#elif OPENFOAMESI
        slip0_.reset(new scalarField(nSlices, 0.0));
#endif
    }

    if (!slip_.valid())
    {
#ifdef OPENFOAMFOUNDATION
        slip_.set(new scalarField(nSlices, 0.0));
#elif OPENFOAMESI
        slip_.reset(new scalarField(nSlices, 0.0));
#endif
    }    
    
    // Update slip0_ once per time-step (used by friction iteration)
    if (timeCheck_ != mesh_.time().timeIndex())
    {
        slip0_() = slip_();
        timeCheck_ = mesh_.time().timeIndex();
    }

    // Material properties
    const volScalarField& mu = mat_.mu();
    const volScalarField& lambda = mat_.lambda();

    // Prepare 2*mu + lambda
    volScalarField twoMuLambda("twoMuLambda", 2*mu + lambda);  
    
    // Additional strains defined by constitutive mechanical laws
    const volSymmTensorField lawsAdditionalStrain(getLawsAdditionalStrain());

    // lambda*(epsElX+epsElY)
    volScalarField lambdaEpsElXEpsElY
    (
        "lambdaEpsElXEpsElY",
        lambda*
        (
            epsilon.component(symmTensor::XX)
          + epsilon.component(symmTensor::YY)
          - additionalStrain_.component(symmTensor::XX)
          - additionalStrain_.component(symmTensor::YY)
          - lawsAdditionalStrain.component(symmTensor::XX)
          - lawsAdditionalStrain.component(symmTensor::YY)
        )
    );

    // (2mu+lambda)*epsAddZ
    volScalarField twoMuLambdaEpsAddZ
    (
        "twoMuLambdaEpsAddZ",
        twoMuLambda*
        (
            additionalStrain_.component(symmTensor::ZZ)
          + lawsAdditionalStrain.component(symmTensor::ZZ)
        )
    );

    // Add thermal expansion contribution
    if (thermalExpansion_)
    {
        twoMuLambdaEpsAddZ += twoMuLambda*epsTh_.component(symmTensor::ZZ);

        lambdaEpsElXEpsElY -= lambda*
        (
            epsTh_.component(symmTensor::XX)
          + epsTh_.component(symmTensor::YY)
        );
    }

    // Old time epsZ (for incremental solvers)
    volScalarField epsilonOldZ
    (
        "epsilonOldZ",
        epsilon.oldTime().component(symmTensor::ZZ)
    );

    // Damage fields
    if (mesh_.foundObject<volTensorField>("damageTensor"))
    {
        const volScalarField damage =
            mesh_.lookupObject<volTensorField>("damageTensor").component(tensor::ZZ);

        twoMuLambda *= (1 - damage);
        lambdaEpsElXEpsElY *= (1 - damage);
        twoMuLambdaEpsAddZ *= (1 - damage);
    }
    else if (mesh_.foundObject<volScalarField>("damage"))
    {
        const volScalarField& damage =
            mesh_.lookupObject<volScalarField>("damage");

        twoMuLambda *= (1 - damage);
        lambdaEpsElXEpsElY *= (1 - damage);
        twoMuLambdaEpsAddZ *= (1 - damage);
    }  

    // Convert 3D fields to slice-average
    const scalarField& twoMuLambda_1D(mapper_->sliceAverage(twoMuLambda));
    const scalarField& lambdaEpsElXEpsElY_1D(mapper_->sliceAverage(lambdaEpsElXEpsElY));
    const scalarField& twoMuLambdaEpsAddZ_1D(mapper_->sliceAverage(twoMuLambdaEpsAddZ));
    const scalarField& epsilonOldZ_1D(mapper_->sliceAverage(epsilonOldZ));

    // Mesh volumes
    const scalarField& V(mesh_.V());

    // Slice type + volumes
    labelList isFuel(nSlices, 0);
    scalarField volumePerSlice(nSlices, 0.0);

    const PtrList<labelList>& sliceAddrList(mapper_->sliceAddrList());

    forAll(sliceAddrList, sliceID)
    {
        const labelList& addr(sliceAddrList[sliceID]);

        isFuel[sliceID] = mapper_->isFuel()[sliceID];

        scalarField Vaddr(V, addr);
        volumePerSlice[sliceID] = gSum(Vaddr);
    }

   // Slice heights (computed from points)
   scalarField sliceHeights(nSlices, 0.0);

   const vectorField& points = mesh_.points();
   const labelListList& cellPoints = mesh_.cellPoints();

   const scalarField z = (points & pinDirection);
   const scalarField r = mag(points - (z*pinDirection));

   scalar fuelTopMaxR(0);
   scalar fuelTopMinR(0);
   scalar cladTopMaxR(0);
   scalar cladTopMinR(0);

   forAll(sliceAddrList, sliceID)
   {
       const labelList& addr(sliceAddrList[sliceID]);

       scalar zMax = -GREAT;
       scalar zMin =  GREAT;

       scalar rMax = -GREAT;
       scalar rMin =  GREAT;

       forAll(addr, addrI)
       {
           const label cellI = addr[addrI];
           const labelList& ptIds = cellPoints[cellI];

           forAll(ptIds, i)
           {
               const label ptI = ptIds[i];
               zMin = min(zMin, z[ptI]);
               zMax = max(zMax, z[ptI]);

               rMax = max(rMax, r[ptI]);
               rMin = min(rMin, r[ptI]);
           }
       }

       reduce(zMin, minOp<scalar>());
       reduce(zMax, maxOp<scalar>());

       reduce(rMin, minOp<scalar>());
       reduce(rMax, maxOp<scalar>());

       if (isFuel[sliceID])
       {
           fuelTopMaxR = rMax;
           fuelTopMinR = rMin;
       }
       else
       {
           cladTopMaxR = rMax;
           cladTopMinR = rMin;
       }

       sliceHeights[sliceID] = zMax - zMin;
   }

    // ---------------------------------------------------------------------
    // Friction prep (gap width, interface pressure, mu, area, old displacements)
    //
    // We cannot just slice-average: we need consistency between opposite sides
    // (fuel/clad) across the coupled interface. Hence AMI mapping on the patch.
    // ---------------------------------------------------------------------
    scalarField gapWidth_1D(nSlices, 0.0);
    scalarField muOwn_1D(nSlices, 0.0);
    scalarField muNbr_1D(nSlices, 0.0);
    scalarField interfacePressure_1D(nSlices, 0.0);
    scalarField magSf_1D(nSlices, 0.0);
    scalarField dispOldOwn_1D(nSlices, 0.0);
    scalarField dispOldNbr_1D(nSlices, 0.0);

    // Will be recomputed inside the force iteration (depends on current D)
    scalarField dispOwn_1D(nSlices, 0.0);
    scalarField dispNbr_1D(nSlices, 0.0);

    // Fields used for friction mapping
    const volScalarField& gapWidth = mesh_.lookupObject<volScalarField>("gapWidth");
    const volScalarField& interfacePressure = mesh_.lookupObject<volScalarField>("interfaceP");
    const volVectorField& dispOld = mesh_.lookupObject<volVectorField>("D").oldTime();
    const surfaceScalarField& magSf = mesh_.magSf();

    forAll(sliceAddrList, sliceID)
    {
        const labelList& sliceAddr = sliceAddrList[sliceID];
        const List<cell> sliceCells(mesh_.cells(), sliceAddr);

        forAll(sliceCells, j)
        {
            const cell& c = sliceCells[j];

            forAll(c, faceI)
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                if (patchID < 0) continue;

                if (!isType<regionCoupledOFFBEATFvPatch>(mesh_.boundary()[patchID]))
                {
                    continue;
                }

                const regionCoupledOFFBEATFvPatch& rcp =
                    refCast<const regionCoupledOFFBEATFvPatch>(mesh_.boundary()[patchID]);

                const label nbrPatchID = rcp.nbrPatchID();

                const scalarField& gapWidthP = gapWidth.boundaryField()[patchID];
                const scalarField& interfacePressureP = interfacePressure.boundaryField()[patchID];
                const scalarField& muOwnP = mu.boundaryField()[patchID];
                const scalarField muNbrP =
                    rcp.regionCoupledPatch().interpolate(mu.boundaryField()[nbrPatchID]);

                const vectorField& dispOldP = dispOld.boundaryField()[patchID];
                const vectorField dispOldNbrP =
                    rcp.regionCoupledPatch().interpolate(dispOld.boundaryField()[nbrPatchID]);

                const scalarField& magSfP = magSf.boundaryField()[patchID];

                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                if (faceID < 0) continue;

                // Contact / overlap convention: gapWidth < 0 => in contact
                if (gapWidthP[faceID] < 0.0)
                {
                    const scalar w = magSfP[faceID];

                    gapWidth_1D[sliceID]         += gapWidthP[faceID]*w;
                    interfacePressure_1D[sliceID] += interfacePressureP[faceID]*w;

                    muOwn_1D[sliceID]            += muOwnP[faceID]*w;
                    muNbr_1D[sliceID]            += muNbrP[faceID]*w;

                    dispOldOwn_1D[sliceID]       += (dispOldP[faceID] & pinDirection)*w;
                    dispOldNbr_1D[sliceID]       += (dispOldNbrP[faceID] & pinDirection)*w;

                    magSf_1D[sliceID]            += w;
                }
            }
        }
    }

    reduce(gapWidth_1D, sumOp<scalarField>());
    reduce(interfacePressure_1D, sumOp<scalarField>());
    reduce(muOwn_1D, sumOp<scalarField>());
    reduce(muNbr_1D, sumOp<scalarField>());
    reduce(dispOldOwn_1D, sumOp<scalarField>());
    reduce(dispOldNbr_1D, sumOp<scalarField>());
    reduce(magSf_1D, sumOp<scalarField>());

    forAll(magSf_1D, sliceID)
    {
        const scalar denom = max(magSf_1D[sliceID], SMALL);
        gapWidth_1D[sliceID]         /= denom;
        interfacePressure_1D[sliceID] /= denom;
        muOwn_1D[sliceID]            /= denom;
        muNbr_1D[sliceID]            /= denom;
        dispOldOwn_1D[sliceID]       /= denom;
        dispOldNbr_1D[sliceID]       /= denom;
    }

    // =========================================================
    // Spring/pressure/friction iteration (axial force balance)
    // =========================================================
    scalar forceAverage(0);
    scalar forceAveragePrev(GREAT);

    scalarField epsilonZ_1D(nSlices, 0.0);

    do
    {
        // Save previous forces for relaxation
        const scalarField frictionForcePrev(frictionForce_());
        const scalarField springForcePrev(springForce_());
        const scalarField pressureForcePrev(pressureForce_());

        forceAveragePrev =
            gSum(mag(springForce_() + pressureForce_() + frictionForce_())*volumePerSlice)
          / gSum(volumePerSlice);

        // Slice-dependent axial strain
        forAll(springForce_(), sliceID)
        {
            const scalar sliceA =
                volumePerSlice[sliceID]/(angularFraction*sliceHeights[sliceID]);

            epsilonZ_1D[sliceID] =
            (
                springForce_()[sliceID]/sliceA
              + pressureForce_()[sliceID]/sliceA
              + frictionForce_()[sliceID]/sliceA
              - lambdaEpsElXEpsElY_1D[sliceID]
              + twoMuLambdaEpsAddZ_1D[sliceID]
            )/(max(twoMuLambda_1D[sliceID], SMALL));
        }

        // Redistribute axial strain to cells (and patch faces)
        scalar totalAxialDispFuel(0);
        scalar totalAxialDispClad(0);
        scalar totalAxialDispFuelOld(0);
        scalar totalAxialDispCladOld(0);

        forAll(sliceAddrList, sliceID)
        {
            const labelList& addr(sliceAddrList[sliceID]);

            scalar totalAxialDisp(0);
            scalar totalAxialDispOld(0);

            if (isFuel[sliceID])
            {
                totalAxialDispFuel += epsilonZ_1D[sliceID]*sliceHeights[sliceID];
                totalAxialDispFuelOld += epsilonOldZ_1D[sliceID]*sliceHeights[sliceID];

                totalAxialDisp = totalAxialDispFuel;
                totalAxialDispOld = totalAxialDispFuelOld;
            }
            else
            {
                totalAxialDispClad += epsilonZ_1D[sliceID]*sliceHeights[sliceID];
                totalAxialDispCladOld += epsilonOldZ_1D[sliceID]*sliceHeights[sliceID];

                totalAxialDisp = totalAxialDispClad;
                totalAxialDispOld = totalAxialDispCladOld;
            }

            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];

                epsilon[cellI].zz() = epsilonZ_1D[sliceID];

                if (mesh_.foundObject<fvMesh>("referenceMesh"))
                {
                    D[cellI].z() = totalAxialDisp - totalAxialDispOld;
                }
                else
                {
                    D[cellI].z() = totalAxialDisp;
                }

                // Correct z-component on boundary faces using cell value.
                // (Cannot use correctBoundaryConditions: would alter x/y components.)
                const cell& c = mesh_.cells()[cellI];

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    if (patchID < 0) continue;

                    const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                    if (faceID < 0) continue;

                    if
                    (
                        !isType<emptyFvPatchSymmTensorField>(epsilon.boundaryField()[patchID])
                    )
                    {
                        epsilon.boundaryFieldRef()[patchID][faceID].zz() =
                            epsilon[cellI].zz();

                        D.boundaryFieldRef()[patchID][faceID].z() =
                            D[cellI].z();
                    }
                }
            }
        }

        // Build AMI-consistent own/nbr axial displacement for friction
        dispOwn_1D = 0.0;
        dispNbr_1D = 0.0;

        forAll(sliceAddrList, sliceID)
        {
            const labelList& sliceAddr = sliceAddrList[sliceID];
            const List<cell> sliceCells(mesh_.cells(), sliceAddr);

            forAll(sliceCells, j)
            {
                const cell& c = sliceCells[j];

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
                    if (patchID < 0) continue;

                    if (!isType<regionCoupledOFFBEATFvPatch>(mesh_.boundary()[patchID]))
                    {
                        continue;
                    }

                    const regionCoupledOFFBEATFvPatch& rcp =
                        refCast<const regionCoupledOFFBEATFvPatch>(mesh_.boundary()[patchID]);

                    const label nbrPatchID = rcp.nbrPatchID();

                    const vectorField& dispOwnP = D.boundaryField()[patchID];
                    const vectorField dispNbrP =
                        rcp.regionCoupledPatch().interpolate(D.boundaryField()[nbrPatchID]);

                    const scalarField& gapWidthP = gapWidth.boundaryField()[patchID];
                    const scalarField& magSfP = magSf.boundaryField()[patchID];

                    const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                    if (faceID < 0) continue;

                    if (gapWidthP[faceID] < 0.0)
                    {
                        const scalar w = magSfP[faceID];

                        dispOwn_1D[sliceID] += (dispOwnP[faceID] & pinDirection)*w;
                        dispNbr_1D[sliceID] += (dispNbrP[faceID] & pinDirection)*w;
                    }
                }
            }
        }

        reduce(dispOwn_1D, sumOp<scalarField>());
        reduce(dispNbr_1D, sumOp<scalarField>());

        forAll(magSf_1D, sliceID)
        {
            const scalar denom = max(magSf_1D[sliceID], SMALL);
            dispOwn_1D[sliceID] /= denom;
            dispNbr_1D[sliceID] /= denom;
        }

        // Compute axial friction force per slice (penalty + Coulomb cap)
        forAll(sliceAddrList, sliceID)
        {
            if (gapWidth_1D[sliceID] < 0.0)
            {
                const scalar sliceA = 2*pi*fuelTopMaxR*sliceHeights[sliceID];

                const scalar muEff = min(muOwn_1D[sliceID], muNbr_1D[sliceID]);

                const scalar stiffness =
                    penaltyScaleFactFriction_*muEff;

                const scalar frictionMax =
                    frictionCoeff_*interfacePressure_1D[sliceID];

                // Relative axial slip increment (own - nbr), accumulated over time-step
                slip_()[sliceID] =
                    slip0_()[sliceID]
                  + (dispOwn_1D[sliceID] - dispOldOwn_1D[sliceID])
                  - (dispNbr_1D[sliceID] - dispOldNbr_1D[sliceID]);

                const scalar frictionTrial = stiffness*slip_()[sliceID] / sliceHeights[sliceID];

                // Optional scaling for the fitting/transition (kept as you intended)
                const scalar fittingPenaltyCoeff =
                    min
                    (
                        scalar(1.0),
                        mag(interfacePressure_1D[sliceID]/frictionSaturationPressure_)
                    );

                const scalar frictionCapped =
                    min(mag(frictionMax), mag(frictionTrial)*fittingPenaltyCoeff);

                const scalar sgn = -frictionTrial/(mag(frictionTrial) + SMALL);

                frictionForce_()[sliceID] = sgn*frictionCapped*sliceA;
            }
            else
            {
                frictionForce_()[sliceID] = 0.0;
                slip_()[sliceID] = 0.0;
            }
        }

        // Update spring and pressure force
        scalar pFluid(0);

#ifdef OPENFOAMESIORFOUNDATION
        if (fluidPressureList_.valid())
        {
            pFluid = fluidPressureList_->value
            (
                mesh_.time().timeToUserTime(mesh_.time().value())
            );
        }
#endif

        // Total spring elongation (fuel relative to clad)
        const scalar springDisp =
            max(totalAxialDispFuel - totalAxialDispClad + springPreCompression_, 0.0);

        forAll(sliceAddrList, sliceID)
        {
            if (isFuel[sliceID])
            {
                springForce_()[sliceID] = -springModulus_*springDisp;

                const scalar rMax = fuelTopMaxR;
                const scalar rMin = fuelTopMinR;

                pressureForce_()[sliceID] =
                    -gapGas.p()*pi*(pow(rMax, 2.0) - pow(rMin, 2.0));
            }
            else
            {
                springForce_()[sliceID] = springModulus_*springDisp;

                const scalar rMax = cladTopMaxR;
                const scalar rMin = cladTopMinR;

                pressureForce_()[sliceID] =
                    gapGas.p()*pi*(pow(rMin, 2.0)) - pFluid*pi*(pow(rMax, 2.0));
            }
        }

        // Relax forces (kept your relax_ usage as-is)
        frictionForce_() =
            frictionForce_()*relax_ + frictionForcePrev*(1 - relax_);

        springForce_() =
            springForce_()*relax_ + springForcePrev*(1 - relax_);

        pressureForce_() =
            pressureForce_()*relax_ + pressureForcePrev*(1 - relax_);

        forceAverage =
            gSum(mag(springForce_() + pressureForce_() + frictionForce_())*volumePerSlice)
          / gSum(volumePerSlice);
    }
    while
    (
        mag(forceAveragePrev - forceAverage)
      > precision_*mag(forceAveragePrev)
    );
}


Foam::scalar Foam::rheologyByMaterial::nextDeltaT()
{
    scalar deltaT(GREAT);

    const labelListList& matAddrList(mat_.matAddrList());

    forAll(mat_.materialsList(), i)
    {
        const labelList& addr(matAddrList[i]);
        deltaT = 
        min(deltaT, lawsList_[i].nextDeltaT(addr));
    }

    return deltaT;
}


Foam::scalar Foam::rheologyByMaterial::nextMacroTime()
{
    scalar macroTime(GREAT);

    forAll(mat_.materialsList(), i)
    {
        macroTime = 
        min(macroTime,lawsList_[i].nextMacroTime());
    }

    return macroTime;
} 


Foam::scalar Foam::rheologyByMaterial::lastMacroTime()
{
    scalar lastMacroTime(GREAT);   

    forAll(mat_.materialsList(), i)
    {
        lastMacroTime = 
        min(lastMacroTime, lawsList_[i].lastMacroTime());
    }

    return lastMacroTime;
} 

// ************************************************************************* //


