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
    const materials& mat,
    const dictionary& rheologyDict
)
:
    rheology(mesh, mat, rheologyDict),
    precision_(0),
    springModulus_(0),
    springPreCompression_(0),
    mapper_(nullptr),
    springForce_(),
    pressureForce_()
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
    // component of the stress and it updates the elastic strain component
    forAll(mat_.materialsList(), i)
    {        
        const labelList& addr(matAddrList[i]);
        lawsList_[i].correct(sigmaHydExplicit, sigmaDev_, epsEl_, addr);
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
    // Apply damage if applicable.
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
    if(mesh_.foundObject<sliceMapper>("sliceMapper"))
    {
        if(mapper_ == nullptr)
        {
            mapper_ = &mesh_.lookupObject<sliceMapper>("sliceMapper");
        }

        // Include the gap gas pressure
        const gapGasModel& gapGas
        = mesh_.lookupObject<gapGasModel>("gapGas");

        // Pin direction assumed to be the z-axis
        vector pinDirection_(0, 0, 1);

        const globalOptions& globalOpt
        (mesh_.lookupObject<globalOptions>("globalOptions"));

        // Obtain model angular fraction to calculate the slice quantities for the 
        // corresponding 360 degrees 3-D rod
        scalar angularFraction = globalOpt.angularFraction();

        // Get number of slices. Slices are equivalent to cellZones.
        int nSlices(mapper_->nSlices());

        if(nSlices<=0)
        {
            FatalErrorInFunction()
            << "When the option \'modifiedPlaneStrain\' is activated, the number "
            << "of slices cannot be less than 1, "<< nl << "as the slices are needed "
            << "to perform the axial force balance." << nl
            << "Check that a sliceMapper different than \"none\" is used "
            << "or check the slice definition."
            << abort(FatalError);
        }

        if(!springForce_.valid())
        {
#ifdef OPENFOAMFOUNDATION            
            springForce_.set(new scalarField(nSlices, 0.0));
#elif OPENFOAMESI            
            springForce_.reset(new scalarField(nSlices, 0.0));
#endif            
        }

        if(!pressureForce_.valid())
        {
#ifdef OPENFOAMFOUNDATION            
            pressureForce_.set(new scalarField(nSlices, 0.0));
#elif OPENFOAMESI            
            pressureForce_.reset(new scalarField(nSlices, 0.0));
#endif   
        }

        // Reference to material properties
        const volScalarField& mu(mat_.mu());
        const volScalarField& lambda(mat_.lambda());

        // Assume that:
        //    sigmaZ =
        //   (2mu+lambda)*epsZ - (2mu+lambda)*epsZadd + lambda*(epsElx+epsElY)
        //  
        //    integral(sigmaZ*dA) = -k_spring*D_top -p_gas*Atop
        //
        // Substituting sigmaZ in the second equaition (and assuming epsZ costant
        // in each slice):
        //    epsZ_slice*int((2mu+lambda)*dA) =
        //    int((2mu+lambda)*epsZadd*dA) - int(lambda*(epsElx+epsElY)*dA) - 
        //    k_spring*D_top -p_gas*Atop
        //    
        // therefore:
        //    epsZ_slice = 1/(int((2mu+lambda)*dA))*
        //    [int((2mu+lambda)*epsZadd*dA) - int(lambda*(epsElx+epsElY)*dA) 
        //      -k_spring*D_top -p_gas*Atop]
        // 
        // Note that the integral terms are substituded with the volume average
        // times the area, eg:
        //         int((2mu+lambda)*dA) = (2mu+lambda)_avg*A 
        //
        // Therefore:
        //    epsZ_slice = 1/((2mu+lambda)_avg*A)*
        //    [ ((2mu+lambda)*epsZadd)_avg*A - (lambda*(epsElx+epsElY))_avg*A 
        //      -k_spring*D_top -p_gas*Atop]
        //
        // Simplifing and dropping the avg subscript:
        //    epsZ_slice = 1/(2mu+lambda)*
        //    [ (2mu+lambda)*epsZadd - lambda*(epsElx+epsElY) 
        //      -k_spring*D_top/A -p_gas*Atop/A]
        //
        // where the area of the slice A can be expressed as:
        //                A = V/(alpha*H)
        //
        // where V is the slice volume, H is the slice Height and alpha is the 
        // angular fraction (e.g. 1.0 for a 360degree 3-D model)

        // Prepare twoMuLambda
        volScalarField twoMuLambda("twoMuLambda", 2*mu + lambda);

        // Prepare lambda*(epsElX+epsElY) field        

        volScalarField lambdaEpsElXEpsElY
        (
            "lambdaEpsElXEpsElY", 
            lambda*
            (
                epsEl_.component(symmTensor::XX) + 
                epsEl_.component(symmTensor::YY)
            )
        );

        // Prepare (2mu+lambda)*epsZadd

        // Get additional strains defined by the constitutive mechanical laws
        volSymmTensorField lawsAdditionalStrain
        (
            getLawsAdditionalStrain()
        );

        volScalarField epsilonAddZ
        (
            additionalStrain_.component(symmTensor::ZZ) 
            + lawsAdditionalStrain.component(symmTensor::ZZ)
        );

        if(thermalExpansion_)
        {
            epsilonAddZ += epsTh_.component(symmTensor::ZZ);
        }

        volScalarField twoMuLambdaEpsAddZ
        (
            "twoMuLambdaEpsAddZ", 
            twoMuLambda*epsilonAddZ
        );

        // Old time epsZ (for incremental solvers)
        volScalarField epsilonOldZ
        (
            "epsilonOldZ", 
            epsilon.oldTime().component(symmTensor::ZZ)
        );
        
        // Convert 3D fields to slice-average
        const scalarField& twoMuLambda_1D(mapper_->sliceAverage(twoMuLambda));
        const scalarField& lambdaEpsElXEpsElY_1D(mapper_->sliceAverage(lambdaEpsElXEpsElY));
        const scalarField& twoMuLambdaEpsAddZ_1D(mapper_->sliceAverage(twoMuLambdaEpsAddZ));
        const scalarField& epsilonOldZ_1D(mapper_->sliceAverage(epsilonOldZ));

        // Reference to mesh volumes
        const scalarField& V(mesh_.V());

        // Prepare list of bool indicating if cellZone is fuel type or not
        labelList isFuel(nSlices, 0);

        // Prepare slice volumes list
        scalarField volumePerSlice(nSlices, 0.0);

        const PtrList<labelList>& sliceAddrList(mapper_->sliceAddrList());

        forAll(sliceAddrList, i)
        {
            const labelList& addr(sliceAddrList[i]);
            // if(addr.size())
            {
                isFuel[i] = mapper_->isFuel()[i];

                scalarField Vaddr(V, addr);
                volumePerSlice[i] = gSum(Vaddr);
            }
        }    

        // Prepare 1D epsilonZ field
        scalarField epsilonZ_1D(nSlices, 0); 

        // Prepare slice heights.
        // It is assumed that the slices for fuel and cladding are given in order
        // from bottom to top. It is also assumed that the cladding caps are not
        // modeled, i.e. there might be a bottom or top plenum with an additional 
        // portion of the cladding outside the bounds of the fuel column, but the 
        // cladding model is truncated
        scalarField sliceHeights(nSlices, 0.0);

        const vectorField& points = mesh_.points();
        const labelListList& cellPoints = mesh_.cellPoints();
        scalarField z = (points & pinDirection_);
        scalarField r = mag(points - (z*pinDirection_));

        scalar fuelTopMaxR(0);
        scalar fuelTopMinR(0);
        scalar cladTopMaxR(0);
        scalar cladTopMinR(0);

        const PtrList<labelList>& sliceAddr(mapper_->sliceAddrList());
        
        // Loop over slices
        forAll(sliceAddr, sliceID)
        {
            const labelList& addr(mapper_->sliceAddrList()[sliceID]);

            scalar zMax = -GREAT;
            scalar zMin = GREAT;

            scalar rMax = -GREAT;
            scalar rMin = GREAT;
            
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

            if(isFuel[sliceID])
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

        // Prepare spring force quantities
        scalar forceAverage(0);
        scalar forceAveragePrev(GREAT);

        // The axial strains are a function of the spring force, which in turn is 
        // calculated as a function of the axial strain. Iterate until convergence.
        do
        {
            forceAveragePrev = 
            gSum((springForce_()+pressureForce_())*volumePerSlice)/gSum(volumePerSlice);

            // Calculate slice-average axial strain
            forAll(springForce_(), sliceID)
            {
                scalar sliceA = 
                volumePerSlice[sliceID]/(angularFraction*sliceHeights[sliceID]);

                epsilonZ_1D[sliceID] = 
                (
                    springForce_()[sliceID]/sliceA +
                    pressureForce_()[sliceID]/sliceA - 
                    lambdaEpsElXEpsElY_1D[sliceID] + 
                    twoMuLambdaEpsAddZ_1D[sliceID]
                )/(max(twoMuLambda_1D[sliceID], SMALL));
            }

            scalar totalAxialDispFuel(0);
            scalar totalAxialDispClad(0);
            scalar totalAxialDispFuelOld(0);
            scalar totalAxialDispCladOld(0);

            // Redistribute average axial strain in cells and faces
            forAll(mapper_->sliceAddrList(), sliceID)
            {
                const labelList addr(mapper_->sliceAddrList()[sliceID]);

                scalar totalAxialDisp(0);
                scalar totalAxialDispOld(0);

                if(isFuel[sliceID])
                {
                    totalAxialDispFuel += 
                    epsilonZ_1D[sliceID]*sliceHeights[sliceID];

                    totalAxialDispFuelOld += 
                    epsilonOldZ_1D[sliceID]*sliceHeights[sliceID];

                    totalAxialDisp = totalAxialDispFuel;
                    totalAxialDispOld = totalAxialDispFuelOld;
                }
                else
                {
                    totalAxialDispClad += 
                    epsilonZ_1D[sliceID]*sliceHeights[sliceID];

                    totalAxialDispCladOld += 
                    epsilonOldZ_1D[sliceID]*sliceHeights[sliceID];

                    totalAxialDisp = totalAxialDispClad;
                    totalAxialDispOld = totalAxialDispCladOld;
                }

                forAll(addr, addrI)
                {
                    const label& cellI(addr[addrI]);

                    epsilon[cellI].zz() = epsilonZ_1D[sliceID];
                    if(mesh_.foundObject<fvMesh>("referenceMesh"))
                    {
                        D[cellI].z() = totalAxialDisp - totalAxialDispOld;
                    }
                    else
                    {
                        D[cellI].z() = totalAxialDisp;
                    }

                    // Correct z component of the strain on the patch face
                    // if current cell is adjacent to boundary patch.
                    // Cannot use correctBoundaryConditions as it alters also
                    // the strain components along x and y

                    // Reference to list of cell faces
                    const cell& c = mesh_.cells()[cellI];     

                    forAll(c, faceI)
                    {
                        const label patchID = 
                        mesh_.boundaryMesh().whichPatch(c[faceI]);

                        if (patchID > -1)
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                            
                            if
                            (
                                !isType<emptyFvPatchSymmTensorField>
                                (
                                    epsilon.boundaryField()[patchID]
                                )
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
            }                   

            // Update value of spring and pressure force 

            scalar pFluid(0);

#ifdef OPENFOAMESIORFOUNDATION           
            if(fluidPressureList_.valid())
            {
                fluidPressureList_->value
                (
                    mesh_.time().timeToUserTime
                    (
                        mesh_.time().value()
                    )
                );
            }
#endif

            scalar pi = Foam::constant::mathematical::pi;

            // Total spring elongation
            scalar springDisp = 
            max(totalAxialDispFuel - totalAxialDispClad + springPreCompression_, 0.0);

            // Calculate spring force for the fuel column
            forAll(mapper_->sliceAddrList(), sliceID)
            {
                if(isFuel[sliceID])
                {
                    springForce_()[sliceID] = -springModulus_*springDisp;

                    scalar rMax = fuelTopMaxR;
                    scalar rMin = fuelTopMinR;

                    pressureForce_()[sliceID] = 
                    -gapGas.p()*pi*(pow(rMax, 2.0) - pow(rMin, 2.0));
                }
                else
                {
                    springForce_()[sliceID]= springModulus_*springDisp;

                    scalar rMax = cladTopMaxR;
                    scalar rMin = cladTopMinR;

                    pressureForce_()[sliceID] = 
                    gapGas.p()*pi*(pow(rMin, 2.0)) - pFluid*pi*(pow(rMax, 2.0));
                }
            }

            forceAverage = 
            gSum((springForce_()+pressureForce_())*volumePerSlice)/gSum(volumePerSlice);
        }
        while
        ( 
            (mag(forceAveragePrev - forceAverage)) > 
            (precision_*mag(forceAveragePrev))
        );
    }
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
