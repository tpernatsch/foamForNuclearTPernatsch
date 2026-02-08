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

#include "misesPlasticity.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(misesPlasticity, 0);
    addToRunTimeSelectionTable
    (
        constitutiveLaw, 
        misesPlasticity, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::misesPlasticity::misesPlasticity
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    constitutiveLaw(mesh, lawDict), 
    yieldStressModel_(yieldStressModel::New(mesh, lawDict)),
    epsilonP_(createOrLookup<symmTensor>(mesh, "epsilonP")),
    DEpsilonP_(createOrLookup<symmTensor>(mesh, "DEpsilonP")),
    sigmaY_(createOrLookup<scalar>(mesh, "sigmaY", dimPressure)),
    DSigmaY_(createOrLookup<scalar>(mesh, "DSigmaY", dimPressure)),
    epsilonPEq_(createOrLookup<scalar>(mesh, "epsilonPEq")),
    DEpsilonPEq_(createOrLookup<scalar>(mesh, "DEpsilonPEq")),
    DLambda_(createOrLookup<scalar>(mesh, "DLambda")),
    activeYield_(createOrLookup<scalar>(mesh, "activeYield")),
    plasticN_(createOrLookup<symmTensor>(mesh, "plasticN")),
    Hp_(0.0),
    maxAveragePlasticIncrement_(GREAT),
    maxMaximumPlasticIncrement_(GREAT)
{    
    epsilonP_.oldTime();
    DEpsilonP_.oldTime();
    sigmaY_.oldTime();
    DSigmaY_.oldTime();
    epsilonPEq_.oldTime();
    DEpsilonPEq_.oldTime();
    DLambda_.oldTime();
    activeYield_.oldTime();
    plasticN_.oldTime();

    const dictionary& controlDict(mesh_.time().controlDict());

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
        maxAveragePlasticIncrement_ = controlDict.lookupOrDefault<scalar>(
            "maxAveragePlasticIncrement", GREAT);
        maxMaximumPlasticIncrement_ = controlDict.lookupOrDefault<scalar>(
            "maxMaximumPlasticIncrement", GREAT);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::misesPlasticity::~misesPlasticity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::misesPlasticity::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{ 
    //- Calculate sigmaY
    yieldStressModel_->correctYieldStress(addr);

    //- TODO: Instead of looping over cells one might store the addressing and
    //- work with openfoam field syntax.
    const volScalarField& mu_(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& K_(mesh_.lookupObject<volScalarField>("threeK"));

    // Hardening slope
    const volScalarField& Hp = yieldStressModel_->Hp();

    // Take references to the internal fields for efficiency
    symmTensorField& plasticNI = plasticN_.ref();
    scalarField& DSigmaYI = DSigmaY_.ref();
    scalarField& DLambdaI = DLambda_.ref();
    const scalarField& muI = mu_.internalField();

    //- Take references to internal fields for fields that must be changed
    scalarField& DEpsilonPEqI = DEpsilonPEq_.ref();
    symmTensorField& DEpsilonPI = DEpsilonP_.ref();
    symmTensorField& epsilonI = epsilonEl.ref();

    symmTensorField& epsilonPI = epsilonP_.ref();
    scalarField& epsilonPEqI =  epsilonPEq_.ref();
    const symmTensorField& epsilonPOldTimeI = epsilonP_.oldTime().internalField(); 
    const scalarField& epsilonPEqOldTimeI = epsilonPEq_.oldTime().internalField();

    //- Check if the simulation is large strain
    // bool largeStrain = mesh_.foundObject<volTensorField>("F") ?
    // true : false;

    scalar relax = mesh_.relaxField(DEpsilonP_.name()) ? 
    mesh_.fieldRelaxationFactor(DEpsilonP_.name()) : 1.0;

    // Calculate DLambda_ and plasticN_
    forAll(addr, addrI)
    {
        //- Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];

        // Calculate deviatoric trial stress
        const symmTensor sTrial = 2.0*mu_[cellI]*dev(epsilonEl[cellI]);
        const scalar magS = mag(sTrial);
        
        // Calculate the yield function
        const scalar fTrial = magS - sqrt(2.0/3.0)*sigmaY_.internalField()[cellI];

        // Calculate return direction plasticN
        if (magS > SMALL)
        {
            plasticNI[cellI] = sTrial/magS;
        }

        // Calculate DLambda/DEpsilonPEq
        if (fTrial < SMALL)
        {
            // elastic
            DSigmaYI[cellI] = 0.0;
            DLambdaI[cellI] = 0.0;
            plasticNI[cellI] = symmTensor::zero;
        }
        else
        {
            // Plastic modulus is linear
            DLambdaI[cellI] = fTrial/(2*muI[cellI]);

            scalar magHp = mag(Hp[cellI]);

            if (magHp > SMALL)
            {
                DLambdaI[cellI] /= 1.0 + Hp[cellI]/(3*muI[cellI]);

                // Update increment of yield stress
                DSigmaYI[cellI] = DLambdaI[cellI]*Hp[cellI];
            }
        }        

        // Update DEpsilonPEq
        DEpsilonPEqI[cellI] = (sqrt(2.0/3.0)*DLambdaI[cellI]);

        // Update DEpsilonP
        DEpsilonPI[cellI] = dev(relax*DLambdaI[cellI]*plasticNI[cellI] + (1 - relax)*DEpsilonPI[cellI]);

        // Update total plastic strain
        epsilonPI[cellI] = epsilonPOldTimeI[cellI] + DEpsilonPI[cellI];

        // Update elastic strain
        epsilonI[cellI] -= DEpsilonPI[cellI];

        // Update equivalent total plastic strain
        epsilonPEqI[cellI] = epsilonPEqOldTimeI[cellI] + DEpsilonPEqI[cellI];

        // Calculate deviatoric stress
        sigmaDev[cellI] = sTrial - 2*muI[cellI]*DEpsilonPI[cellI];

        // Calculate the hydrostatic pressure
        const scalar trEpsilon = tr(epsilonI[cellI]);

        sigmaHyd[cellI] = 1.0/3.0*K_[cellI]*trEpsilon;

        // epsilonI[cellI] -= DEpsilonPI[cellI];

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            // forAll(patchIDs, registeredID)
            {
                if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);           
    
                    // Take references to the boundary patch fields for efficiency
                    symmTensorField& sigmaDevP = 
                    sigmaDev.boundaryFieldRef()[patchID];

                    scalarField& sigmaHydP = 
                    sigmaHyd.boundaryFieldRef()[patchID];

                    symmTensorField& plasticNP = 
                    plasticN_.boundaryFieldRef()[patchID];

                    scalarField& DSigmaYP = 
                    DSigmaY_.boundaryFieldRef()[patchID];

                    scalarField& DLambdaP = 
                    DLambda_.boundaryFieldRef()[patchID];

                    const scalarField& muP = 
                    mu_.boundaryField()[patchID];

                    const scalarField& sigmaYP = 
                    sigmaY_.boundaryField()[patchID];
                    const scalarField& HpP = 
                    Hp.boundaryField()[patchID];

                    scalarField& DEpsilonPEqP = 
                    DEpsilonPEq_.boundaryFieldRef()[patchID];

                    symmTensorField& DEpsilonPP = 
                    DEpsilonP_.boundaryFieldRef()[patchID];

                    symmTensorField& epsilonP = 
                    epsilonEl.boundaryFieldRef()[patchID];

                    symmTensorField& epsilonPP = 
                    epsilonP_.boundaryFieldRef()[patchID];

                    scalarField& epsilonPEqP = 
                    epsilonPEq_.boundaryFieldRef()[patchID];

                    const scalarField& KP = K_.boundaryField()[patchID];

                    const symmTensorField& epsilonPOldTimeP = 
                    epsilonP_.oldTime().boundaryField()[patchID]; 
                    const scalarField& epsilonPEqOldTimeP = 
                    epsilonPEq_.oldTime().boundaryField()[patchID]; 

                    // Calculate deviatoric trial stress
                    const symmTensor sTrialP = 
                    2.0*muP[faceID]*(dev(epsilonP[faceID]));

                    // Calculate the yield function
                    const scalar fTrialP = 
                    mag(sTrialP) - sqrt(2.0/3.0)*sigmaYP[faceID];

                    // Calculate direction plasticN
                    const scalar magS = mag(sTrialP);
                    if (magS > SMALL)
                    {
                        plasticNP[faceID] = sTrialP/magS;
                    }

                    // Calculate DEpsilonPEq
                    if (fTrialP < SMALL)
                    {
                        // elasticity
                        DSigmaYP[faceID] = 0.0;
                        DLambdaP[faceID] = 0.0;
                        plasticNP[faceID] = symmTensor::zero;
                    }
                    else
                    {
                        // Plastic modulus is linear
                        DLambdaP[faceID] = fTrialP/(2.0*muP[faceID]);

                        scalar magHp = mag(HpP[faceID]);
                        if (magHp > SMALL)
                        {
                            DLambdaP[faceID] /= 1.0 + HpP[faceID]/(3.0*muP[faceID]);

                            // Update increment of yield stress
                            DSigmaYP[faceID] = DLambdaP[faceID]*HpP[faceID];
                        }
                    }

                    // Update DEpsilonPEq
                    DEpsilonPEqP[faceID] = (sqrt(2.0/3.0)*DLambdaP[faceID]);

                    // Update DEpsilonP
                    DEpsilonPP[faceID] = dev(relax*DLambdaP[faceID]*plasticNP[faceID] + (1 - relax)*DEpsilonPP[faceID]);

                    // Update total plastic strain
                    epsilonPP[faceID] = 
                    epsilonPOldTimeP[faceID] + DEpsilonPP[faceID];

                    // Update elastic strain
                    epsilonP[faceID] -=DEpsilonPP[faceID];

                    // Update equivalent total plastic strain
                    epsilonPEqP[faceID] = 
                    epsilonPEqOldTimeP[faceID] + DEpsilonPEqP[faceID];

                    // Calculate deviatoric stress
                    sigmaDevP[faceID] = 
                    sTrialP - 2*muP[faceID]*DEpsilonPP[faceID];

                    // Calculate the hydrostatic pressure
                    const scalar trEpsilonP = tr(epsilonP[faceID]);

                    sigmaHydP[faceID] = 1.0/3.0*KP[faceID]*trEpsilonP;
                }
            }
        }
    }
}


void Foam::misesPlasticity::correctEpsilonEl
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        //- Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];

        epsilonEl[cellI] -= epsilonP_.oldTime()[cellI];
        
        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            {
                if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);     

                    // Take references to the boundary patch fields        
                    symmTensorField& epsilonEl_p = 
                    epsilonEl.boundaryFieldRef()[patchID]; 
                    
                    const symmTensorField& epsilonPOldTime_p = 
                    epsilonP_.oldTime().boundaryField()[patchID]; 

                    epsilonEl_p[faceID] -= epsilonPOldTime_p[faceID];
                }
            }
        }
    }
};


void Foam::misesPlasticity::updateTotalFields
(
    const labelList& addr
)
{
    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        sigmaY_.ref()[cellI] += DSigmaY_.internalField()[cellI];

        if (DEpsilonPEq_.internalField()[cellI] > SMALL)
        {
            activeYield_.ref()[cellI] = 1.0;
            numCellsYielding++;
        }
        else
        {
            activeYield_.ref()[cellI] = 0.0;
        }

        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            if 
            (
                patchID > -1 
                and DEpsilonPEq_.boundaryField()[patchID].size()
                and !activeYield_.boundaryField()[patchID].coupled()
            )
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                sigmaY_.boundaryFieldRef()[patchID][faceID] += 
                DSigmaY_.boundaryFieldRef()[patchID][faceID];

                if (DEpsilonPEq_.boundaryField()[patchID][faceID] > SMALL)
                {
                    activeYield_.boundaryFieldRef()[patchID][faceID] = 1.0;
                }
                else
                {
                    activeYield_.boundaryFieldRef()[patchID][faceID] = 0.0;
                }

            }
        }
    }

    reduce(numCellsYielding, sumOp<int>());

    activeYield_.correctBoundaryConditions();

    const int nTotalCells = returnReduce(mesh_.nCells(), sumOp<int>());

    Info<< "    " << numCellsYielding << " cells ("
        << 100.0*scalar(numCellsYielding)/scalar(nTotalCells)
        << "% of the cells in this material) are actively yielding"
        << nl << endl;
    
}


void Foam::misesPlasticity::correctAdditionalStrain
(
    volSymmTensorField& additionalStrain,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        additionalStrain[cellI] = epsilonP_.internalField()[cellI];
    }
}


Foam::scalar Foam::misesPlasticity::nextDeltaT
(
    const labelList& addr
)
{
    //- Store old value of deltaT
    scalar oldDeltaT = mesh_.time().deltaT().value();

    //- Prepare new deltaT
    scalar deltaT(0);

    //- Prepare average and maximum DepsilonPEq from last time step
    scalar averageDepsilonPEq(0);
    scalar maxDepsilonPEq(0);
    scalar volTot(0);

    const scalarField& DepsilonPEq(DEpsilonPEq_.internalField());
    const scalarField& V(mesh_.V());

    //- Find maximum and average DepsilonCreepEq
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        averageDepsilonPEq += (DepsilonPEq[cellI])*V[cellI];
        volTot += V[cellI];

        maxDepsilonPEq = max(maxDepsilonPEq, DepsilonPEq[cellI]);
    }

    reduce(volTot, sumOp<scalar>());
    reduce(averageDepsilonPEq, sumOp<scalar>());
    reduce(maxDepsilonPEq, maxOp<scalar>());

    averageDepsilonPEq /= volTot;

    //- Calculate average and maximum creep rate from last time step
    scalar averagePlasticStrainRate(averageDepsilonPEq/oldDeltaT);
    scalar maxPlasticStrainRate(maxDepsilonPEq/oldDeltaT);

    deltaT = min
    (
        maxAveragePlasticIncrement_/max(averagePlasticStrainRate, SMALL), 
        maxMaximumPlasticIncrement_/max(maxPlasticStrainRate, SMALL)
    );
    
    Info<< "Maximum deltaT calculated by misesPlasticity model: " 
    << mesh_.time().timeToUserTime(deltaT)
    << " with a maximum plastic strain rate of " << maxPlasticStrainRate 
    << endl;

    return deltaT;
}

// ************************************************************************* //


