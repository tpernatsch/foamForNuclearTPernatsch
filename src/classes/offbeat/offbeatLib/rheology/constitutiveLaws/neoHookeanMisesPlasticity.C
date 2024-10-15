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

#include "neoHookeanMisesPlasticity.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "Tuple2.H"
#include "PstreamReduceOps.H"
#include "transformGeometricField.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(neoHookeanMisesPlasticity, 0);
    addToRunTimeSelectionTable
    (
        constitutiveLaw, 
        neoHookeanMisesPlasticity, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::neoHookeanMisesPlasticity::neoHookeanMisesPlasticity
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    constitutiveLaw(mesh, lawDict), 
    sigmaY_(createOrLookup<scalar>(mesh, "sigmaY", dimPressure)),
    DSigmaY_(createOrLookup<scalar>(mesh, "DSigmaY", dimPressure)),
    epsilonP_(createOrLookup<symmTensor>(mesh, "epsilonP")),
    DEpsilonP_(createOrLookup<symmTensor>(mesh, "DEpsilonP")),
    epsilonPEq_(createOrLookup<scalar>(mesh, "epsilonPEq")),
    DEpsilonPEq_(createOrLookup<scalar>(mesh, "DEpsilonPEq")),
    DLambda_(createOrLookup<scalar>(mesh, "DLambda")),
    bEbarTrial_(createOrLookup(mesh, "bEbarTrial", dimless, symmTensor::I)),
    bEbar_(createOrLookup(mesh, "bEbar", dimless, symmTensor::I)),
    plasticN_(createOrLookup<symmTensor>(mesh, "plasticN")),
    activeYield_(createOrLookup<scalar>(mesh, "activeYield")),
    Hp_(0.0) 
{  
#ifdef OPENFOAMFOUNDATION
    stressPlasticStrainSeries_.set(new Function1s::Table<scalar>
        (
            "plasticStrainVsYieldStress", lawDict
        )
    );
#elif OPENFOAMESI
    stressPlasticStrainSeries_.reset(new Table
        (
            "plasticStrainVsYieldStress", lawDict
        )
    );
#endif     
    epsilonP_.oldTime();
    DEpsilonP_.oldTime();
    sigmaY_.oldTime();
    DSigmaY_.oldTime();
    epsilonPEq_.oldTime();
    DEpsilonPEq_.oldTime();
    DLambda_.oldTime();
    bEbarTrial_.oldTime();
    bEbar_.oldTime();
    plasticN_.oldTime();
    activeYield_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::neoHookeanMisesPlasticity::~neoHookeanMisesPlasticity()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::neoHookeanMisesPlasticity::correct
(
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev, 
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    //- Adjust sigmaY if it is first time (sigmaY = 0.0)
    forAll(addr, addrI)
    {
        //- Cell adrress. Act only on element of addr
        const label cellI = addr[addrI];

        const scalarField& epsilonPEqI(epsilonPEq_.internalField());

        if(sigmaY_.ref()[cellI] == 0.0)
        {
            sigmaY_.ref()[cellI] = 
            stressPlasticStrainSeries_->value(epsilonPEqI[cellI]);
        }

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = 
            mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                scalarField& sigmaYP(sigmaY_.boundaryFieldRef()[patchID]);
                
                const scalarField& epsilonPEqP
                (
                    epsilonPEq_.boundaryField()[patchID]
                );

                if(sigmaYP[faceID] == 0.0)
                {
                    sigmaYP[faceID] = 
                    stressPlasticStrainSeries_->value(epsilonPEqP[faceID]);
                }
            }
        }
    }

    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& threeK(mesh_.lookupObject<volScalarField>("threeK"));
    const volScalarField& J(mesh_.lookupObject<volScalarField>("J"));
    const volScalarField& JOld(mesh_.lookupObject<volScalarField>("J").oldTime());

    // TODO (EB) : I created the field relF in the largeStrainTot class as it 
    //             needs to be used here, even if not used anywhere else.
    const volTensorField& relF(mesh_.lookupObject<volTensorField>("relF"));

    // Take references to the internal fields for efficiency
    symmTensorField&     plasticNI = plasticN_.ref();
    symmTensorField&   bEbarTrialI = bEbarTrial_.ref();
    symmTensorField&        bEbarI = bEbar_.ref();
    const symmTensorField& bEbarOldTimeI = bEbar_.oldTime().internalField();
    scalarField&           sigmaYI = sigmaY_.ref();
    scalarField&          DSigmaYI = DSigmaY_.ref();
    scalarField&          DLambdaI = DLambda_.ref();
    const scalarField&         muI = mu.internalField();

    symmTensorField&     epsilonPI = epsilonP_.ref();
    scalarField&      DEpsilonPEqI = DEpsilonPEq_.ref();
    symmTensorField&    DEpsilonPI = DEpsilonP_.ref();
    scalarField&       epsilonPEqI =  epsilonPEq_.ref();
    const symmTensorField& epsilonPOldTimeI = 
        epsilonP_.oldTime().internalField(); 
    const scalarField&   epsilonPEqOldTimeI = 
        epsilonPEq_.oldTime().internalField();

    forAll(addr, addrI)
    {

        //- Cell address. Act only on element of addr
        const label cellI = addr[addrI];

        const scalar relJ(J[cellI]/JOld[cellI]);

        // Calculate the relative deformation gradient with the volumetric term
        // removed
        const tensor relFbar = pow(relJ, -1.0/3.0)*relF[cellI];

        // Update bE trial
        bEbarTrialI[cellI] = transform(relFbar, bEbarOldTimeI[cellI]);

        // Calculate trial deviatoric stress
        const symmTensor sTrial(muI[cellI]*dev(bEbarTrialI[cellI]));

        const scalar Ibar(tr(bEbarTrialI[cellI])/3.0);
        const scalar muBar(Ibar*muI[cellI]);

        // Check for plastic loading and calculate increment of plastic 
        // equivalent strain i.e. the plastic multiplier

        // Trial yield function
        // sigmaY is the Cauchy yield stress so we scale it by J
        const scalar fTrial = mag(sTrial) - sqrt(2.0/3.0)*J[cellI]*sigmaYI[cellI];

        // Magnitude of hardening slope
        const scalar magHp = mag(Hp_);

        // Calculate return direction plasticN
        const scalar magS = mag(sTrial);
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
        }
        else
        {
            // Plastic modulus is linear
            DLambdaI[cellI] = fTrial/(2*muBar);

            if (magHp > SMALL)
            {
                DLambdaI[cellI] /= 1.0 + Hp_/(3*muBar);

                // Update increment of yield stress
                DSigmaYI[cellI] = DLambdaI[cellI]*Hp_;
            }
        }

        // Update DEpsilonP and DEpsilonPEq
        DEpsilonPEqI[cellI] = sqrt(2.0/3.0)*DLambdaI[cellI];
        DEpsilonPI[cellI] = Ibar*DLambdaI[cellI]*plasticNI[cellI];

        // Update epsilonP
        epsilonPI[cellI]   = epsilonPOldTimeI[cellI]   + DEpsilonPI[cellI];
        epsilonPEqI[cellI] = epsilonPEqOldTimeI[cellI] + DEpsilonPEqI[cellI];

        // Calculate deviatoric stress
        sigmaDev[cellI] = 
        sTrial - 2*mu[cellI]*DEpsilonPI[cellI];

        // Update bEbar
        bEbarI[cellI] = (sigmaDev[cellI]/muI[cellI]) + Ibar*I;

        // Update hydrostatic stress 
        sigmaHyd[cellI] = 
        0.5*threeK[cellI]/3.0*(pow(J[cellI], 2.0) - 1.0);

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]); 

                // Take references to the boundary patch fields for efficiency
                const tensorField& relFP = relF.boundaryField()[patchID];

                const scalarField& JP = J.boundaryField()[patchID];
                const scalarField& JOldP = JOld.boundaryField()[patchID];
                const scalar relJP(JP[faceID]/JOldP[faceID]);

                symmTensorField& plasticNP = 
                plasticN_.boundaryFieldRef()[patchID];
                symmTensorField& bEbarTrialP = 
                bEbarTrial_.boundaryFieldRef()[patchID];
                symmTensorField& bEbarP = 
                bEbar_.boundaryFieldRef()[patchID];
                const symmTensorField& bEbarOldTimeP = 
                bEbar_.oldTime().boundaryField()[patchID];
                scalarField& sigmaYP = sigmaY_.boundaryFieldRef()[patchID];
                scalarField& DSigmaYP = DSigmaY_.boundaryFieldRef()[patchID];
                scalarField& DLambdaP = DLambda_.boundaryFieldRef()[patchID];
                const scalarField& muP = mu.boundaryField()[patchID];
                const scalarField& KP = threeK.boundaryField()[patchID];

                symmTensorField& epsilonPP = 
                epsilonP_.boundaryFieldRef()[patchID];
                scalarField& DEpsilonPEqP = 
                DEpsilonPEq_.boundaryFieldRef()[patchID];
                symmTensorField& DEpsilonPP = 
                DEpsilonP_.boundaryFieldRef()[patchID];
                scalarField& epsilonPEqP =  
                epsilonPEq_.boundaryFieldRef()[patchID];

                symmTensorField& sigmaDevP = sigmaDev.boundaryFieldRef()[patchID];
                scalarField& sigmaHydP = sigmaHyd.boundaryFieldRef()[patchID];

                const symmTensorField& epsilonPOldTimeP = 
                    epsilonP_.oldTime().boundaryFieldRef()[patchID];
                const scalarField&   epsilonPEqOldTimeP = 
                    epsilonPEq_.oldTime().boundaryFieldRef()[patchID];

                // Calculate the relative deformation gradient with the volumetric term
                // removed
                const tensor relFbarP = pow(relJP, -1.0/3.0)*relFP[faceID];

                // Update bE trial
                bEbarTrialP[faceID] = transform(relFbarP, bEbarOldTimeP[faceID]);

                // Calculate trial deviatoric stress
                const symmTensor sTrialP(muP[faceID]*dev(bEbarTrialP[faceID]));

                const scalar IbarP(tr(bEbarTrialP[faceID])/3.0);
                const scalar muBarP(IbarP*muP[faceID]);

                // Check for plastic loading and calculate increment of plastic 
                // equivalent strain i.e. the plastic multiplier

                // Trial yield function
                // sigmaY is the Cauchy yield stress so we scale it by J
                const scalar fTrialP(mag(sTrialP) - sqrt(2.0/3.0)*JP[faceID]*sigmaYP[faceID]);

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
                }
                else
                {
                    // yielding
                    // Plastic modulus is linear
                    DLambdaP[faceID] = fTrialP/(2.0*muBarP);

                    if (magHp > SMALL)
                    {
                        DLambdaP[faceID] /= 1.0 + Hp_/(3.0*muBarP);

                        // Update increment of yield stress
                        DSigmaYP[faceID] = DLambdaP[faceID]*Hp_;
                    }
                }

                // Update DEpsilonP and DEpsilonPEq
                DEpsilonPEqP[faceID] = sqrt(2.0/3.0)*DLambdaP[faceID];
                DEpsilonPP[faceID] = IbarP*DLambdaP[faceID]*plasticNP[faceID];

                // Update epsilonP
                epsilonPP[faceID]   = epsilonPOldTimeP[faceID]   + DEpsilonPP[faceID];
                epsilonPEqP[faceID] = epsilonPEqOldTimeP[faceID] + DEpsilonPEqP[faceID];
                
                // Calculate deviatoric stress
                sigmaDevP[faceID] = 
                sTrialP - 2*muP[faceID]*DEpsilonPP[faceID];

                // Update bEbar
                bEbarP[faceID] = (sigmaDevP[faceID]/muP[faceID]) + IbarP*I;

                // Update hydrostatic stress 
                sigmaHydP[faceID] = 
                0.5*KP[faceID]/3.0*(pow(JP[faceID], 2.0) - 1.0);
            }
        }
    }
}



void Foam::neoHookeanMisesPlasticity::updateStress
(
    volSymmTensorField& sigma, 
    volScalarField& sigmaHyd, 
    volSymmTensorField& sigmaDev,
    const labelList& addr
)
{
    const volScalarField& J(mesh_.lookupObject<volScalarField>("J"));

    damageModel_->correct(addr);
    const volScalarField& d(damageModel_->d());
    bool applyDamage(!isNull(d));

    forAll(addr, addrI)
    {

        //- Cell address. Act only on element of addr
        const label cellI = addr[addrI];

        // Apply damage only when applicable
        if(applyDamage)
        {
            sigmaHyd[cellI] *= (1 - d[cellI]);
            sigmaDev[cellI] *= (1 - d[cellI]);
        }

        // Update the Cauchy stress
        sigma[cellI] = 1.0/J[cellI]*(sigmaHyd[cellI]*I + sigmaDev[cellI]);

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and J.boundaryField()[patchID].size())
            {
                const label faceID = 
                mesh_.boundaryMesh()[patchID].whichFace(c[faceI]); 

                // Take references to the boundary patch fields for efficiency
                const scalarField& JP = J.boundaryField()[patchID];

                symmTensorField& sigmaDevP = 
                sigmaDev.boundaryFieldRef()[patchID];
                scalarField& sigmaHydP = 
                sigmaHyd.boundaryFieldRef()[patchID];
                symmTensorField& sigmaP = 
                sigma.boundaryFieldRef()[patchID];

                // Apply damage only when applicable
                if(applyDamage)
                {
                    const scalarField& dP = 
                    d.boundaryField()[patchID];

                    sigmaHydP[faceID] *= (1 - dP[faceID]);
                    sigmaDevP[faceID] *= (1 - dP[faceID]);
                }

                // Update the Cauchy stress
                sigmaP[faceID] = 1.0/JP[faceID]*
                (sigmaHydP[faceID]*I + sigmaDevP[faceID]);
            }
        }
    }
}

void Foam::neoHookeanMisesPlasticity::updateTotalFields
(
    const labelList& addr
)
{
    // Info<< nl << "Updating total accumulated fields" << endl;
    // Info<< "    Max DEpsilonPEq is " << gMax(DEpsilonPEq_) << endl;

    // Count cells actively yielding
    int numCellsYielding = 0;

    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        sigmaY_[cellI] += DSigmaY_[cellI];

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


void Foam::neoHookeanMisesPlasticity::correctAdditionalStrain
(
    volSymmTensorField& additionalStrain,
    const labelList& addr
)
{
    forAll(addr, addrI)
    {
        const label cellI(addr[addrI]);

        additionalStrain[cellI] = epsilonP_[cellI];
    }
}

// ************************************************************************* //
