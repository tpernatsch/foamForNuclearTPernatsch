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

#include "yieldStressFraptran.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(yieldStressFraptran, 0);
    addToRunTimeSelectionTable
    (
        yieldStressModel, 
        yieldStressFraptran, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::yieldStressFraptran::yieldStressFraptran
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    yieldStressModel(mesh, lawDict),
    phiValue_(readScalar(lawDict.lookup("phi"))),
    CWValue_(readScalar(lawDict.lookup("CW"))),    
    Kcoeff_(createOrLookup<scalar>(mesh, "Kcoeff", dimless, 0.0, "zeroGradient")),
    m_(createOrLookup<scalar>(mesh, "m", dimless, 0.0, "zeroGradient")),
    n_(createOrLookup<scalar>(mesh, "n", dimless, 0.0, "zeroGradient"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yieldStressFraptran::~yieldStressFraptran()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::yieldStressFraptran::correctYieldStress
(
    const labelList& addr
)
{
    ///- Elastic modulus and temperature
    const volScalarField& E(mesh_.lookupObject<volScalarField>("E"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    const volSymmTensorField& epsilon(mesh_.lookupObject<volSymmTensorField>("epsilon"));
    const volSymmTensorField& epsilonOld(mesh_.lookupObject<volSymmTensorField>("epsilon").oldTime());
    const volSymmTensorField& sigma(mesh_.lookupObject<volSymmTensorField>("sigma"));

    phiValue_ *= 1e25;
    //- strain rate exponent m
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        scalar T = Temp.internalField()[cellI];

        if(T < 750)
        {
            m_[cellI] = 0.015;
        }
        if(750 <= T && T <= 800)
        {
            m_[cellI] = 7.458e-4*T - 0.544338;
        }   
        if(T > 800)
        {
            m_[cellI] = 3.24124e-4*T - 0.20701;
        }

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            
            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                
                const scalarField& Tp = Temp.boundaryField()[patchID];
                scalarField& mp = m_.boundaryFieldRef()[patchID];

                T = Tp[faceID];
                
                if(T < 750)
                {
                    mp[faceID] = 0.015;
                }
                if(750 <= T && T <= 800)
                {
                    mp[faceID] = 7.458e-4*T - 0.544338;
                }   
                if(T > 800)
                {
                    mp[faceID] = 3.24124e-4*T - 0.20701;
                }
            }
        }

    }

    //- strain hardening exponent n
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        scalar T = Temp.internalField()[cellI];

        if(T < 419.4)
        {
            n_[cellI] = 0.11405;
        }
        else if(419.4 <= T && T <= 1099.0722)
        {
            n_[cellI] = -9.490e-2 + 1.165e-3*T - 1.992e-6*pow(T, 2) + 9.558e-10*pow(T, 3);
        }
        else if(1099.0722 <= T && T <= 1600)
        {
            n_[cellI] = -0.22655119 + 2.5e-4*T;
        }
        else if(T > 1600)
        {
            n_[cellI] = 0.17344880;
        }

        if(phiValue_ < 0.1e25){
            n_[cellI] *= 1.321 + 0.48e-25*phiValue_;
        }
        else if(0.1e25 <= phiValue_ && phiValue_ < 2e25){
            n_[cellI] *= 1.369 + 0.096e-25*phiValue_;
        }
        else if(2e25 <= phiValue_ && phiValue_ < 7.5e25){
            n_[cellI] *= 1.5435 + 0.008727e-25*phiValue_;
        }
        else if(phiValue_ >= 7.5e25){
            n_[cellI] *= 1.608953;
        }
        
        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            
            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                
                const scalarField& Tp = Temp.boundaryField()[patchID];
                scalarField& np = n_.boundaryFieldRef()[patchID];
                
                T = Tp[faceID];

                if(T < 419.4)
                {
                    np[faceID] = 0.11405;
                }
                else if(419.4 <= T && T <= 1099.0722)
                {
                    np[faceID] = -9.490e-2 + 1.165e-3*T - 1.992e-6*pow(T, 2) + 9.558e-10*pow(T, 3);
                }
                else if(1099.0722 <= T && T <= 1600)
                {
                    np[faceID] = -0.22655119 + 2.5e-4*T;
                }
                else if(T > 1600)
                {
                    np[faceID] = 0.17344880;
                }

                if(phiValue_ < 0.1e25){
                    np[faceID] *= 1.321 + 0.48e-25*phiValue_;
                }
                else if(0.1e25 <= phiValue_ && phiValue_ < 2e25){
                    np[faceID] *= 1.369 + 0.096e-25*phiValue_;
                }
                else if(2e25 <= phiValue_ && phiValue_ < 7.5e25){
                    np[faceID] *= 1.5435 + 0.008727e-25*phiValue_;
                }
                else if(phiValue_ >= 7.5e25){
                    np[faceID] *= 1.608953;
                }
            }
        }
    }

    //- Strength coefficient
    const scalar KCW = 0.546 * CWValue_;
    scalar Kphi = 0;
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        scalar T = Temp.internalField()[cellI];

        if(T < 750)
        {
            Kcoeff_[cellI] = 1.17628e9 + 4.54859e5*T - 3.28185e3*pow(T, 2) + 1.72752*pow(T,3);
        }
        else if(750 <= T && T < 1090)
        {
            Kcoeff_[cellI] = 2.522488e6 * exp(2.8500027e6/pow(T, 2));
        }
        else if(1090 <= T && T < 1255)
        {
            Kcoeff_[cellI] = 1.841376039e8 -  1.4345448e5*T;
        }
        else if(1255 <= T && T < 2100)
        {
            Kcoeff_[cellI] = 4.330e7 - 6.685e4*T + 37.579*pow(T, 2) - 7.33e-3*pow(T, 3);
        }

        if(phiValue_ < 0.1e25){
            Kphi = -0.1461 + 1.464e-25*phiValue_ * ((2.25 * exp(-20*CWValue_) * min(1, exp((T-550)/10))) + 1);
            Kcoeff_[cellI] *= 1 + KCW + Kphi;
        }
        else if(phiValue_ <= 0.1e25 && phiValue_ < 2e25){
            Kphi = 2.928e-26*phiValue_;
            Kcoeff_[cellI] *= 1 + KCW + Kphi;
        }
        else if(phiValue_ <= 2e25 && phiValue_ < 12e25){
            Kphi = 0.53236 + 2.6618e-27*phiValue_;
            Kcoeff_[cellI] *= 1 + KCW + Kphi;
        }
        
        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            
            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                
                const scalarField& Tp = Temp.boundaryField()[patchID];
                scalarField& KcoeffP = Kcoeff_.boundaryFieldRef()[patchID];
                
                T = Tp[faceID];

                if(T < 750)
                {
                    KcoeffP[faceID] = 1.17628e9 + 4.54859e5*T - 3.28185e3*pow(T, 2) + 1.72752*pow(T,3);
                }
                else if(750 <= T && T < 1090)
                {
                    KcoeffP[faceID]  = 2.522488e6 * exp(2.8500027e6/pow(T, 2));
                }
                else if(1090 <= T && T < 1255)
                {
                    KcoeffP[faceID]  = 1.841376039e8 -  1.4345448e5*T;
                }
                else if(1255 <= T && T < 2100)
                {
                    KcoeffP[faceID]  = 4.330e7 - 6.685e4*T + 37.579*pow(T, 2) - 7.33e-3*pow(T, 3);
                }

                if(phiValue_ < 0.1e25){
                    Kphi = -0.1461 + 1.464e-25*phiValue_ * ((2.25 * exp(-20*CWValue_) * min(1, exp((T-550)/10))) + 1);
                    KcoeffP[faceID]  *= 1 + KCW + Kphi;
                }
                else if(phiValue_ <= 0.1e25 && phiValue_ < 2e25){
                    Kphi = 2.928e-26*phiValue_;
                    KcoeffP[faceID]  *= 1 + KCW + Kphi;
                }
                else if(phiValue_ <= 2e25 && phiValue_ < 12e25){
                    Kphi = 0.53236 + 2.6618e-27*phiValue_;
                    KcoeffP[faceID]  *= 1 + KCW + Kphi;
                }
            }
        }
    }

    scalar sigmaEq = 0;
    scalar epsilonEq = 0;
    scalar epsilonEqOld = 0;
    scalar ddtEpsilonEq = 0;

    //- yieldStress update
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        const scalar Ei = E.internalField()[cellI];

        sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigma[cellI])));
        epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilon[cellI])));
        epsilonEqOld = sqrt((2.0/3.0)*magSqr(dev(epsilonOld[cellI])));
        ddtEpsilonEq = max(mag(epsilonEq - epsilonEqOld)/mesh_.time().deltaT0().value(), 1e-3); 

        scalar term1 = Kcoeff_[cellI] / pow(Ei, n_[cellI]) * pow(ddtEpsilonEq/1e-3, m_[cellI]);
        // scalar term1 = Kcoeff_[cellI] / pow(Ei, n_[cellI]) * pow(1, m_[cellI]);        
        scalar term2 = 1/(1-n_[cellI]);
        scalar sigmaYOld = Kcoeff_[cellI]*pow(epsilonEqOld, n_[cellI]) * pow(ddtEpsilonEq/1e-3, m_[cellI]);     
        
        sigmaY_[cellI] = max(pow(term1, term2), sigmaYOld) ;
        
        Hp_[cellI] = max(1/((1/n_[cellI]/Kcoeff_[cellI])*pow(sigmaEq/Kcoeff_[cellI],1/n_[cellI] - 1)*pow(1e-3/ddtEpsilonEq, m_[cellI]/n_[cellI]) - 1/Ei), 0.0);
    }     

    //- yieldStress update
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);
            
            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                const scalarField& EiP = E.boundaryField()[patchID];
                const scalarField& TempP = Temp.boundaryField()[patchID];
                const scalarField& KcoeffP = Kcoeff_.boundaryField()[patchID];
                const scalarField& nP = n_.boundaryField()[patchID];
                const scalarField& mP = m_.boundaryField()[patchID];

                const symmTensorField& sigmaP = sigma.boundaryField()[patchID];
                const symmTensorField& epsilonP = epsilon.boundaryField()[patchID];
                const symmTensorField& epsilonOldP = epsilonOld.boundaryField()[patchID];

                scalarField& sigmaYP(sigmaY_.boundaryFieldRef()[patchID]);                                
                scalarField& HpP(Hp_.boundaryFieldRef()[patchID]);        
                
                sigmaEq = sqrt((3.0/2.0)*magSqr(dev(sigmaP[faceID])));
                epsilonEq = sqrt((2.0/3.0)*magSqr(dev(epsilonP[faceID])));
                epsilonEqOld = sqrt((2.0/3.0)*magSqr(dev(epsilonOldP[faceID])));
                ddtEpsilonEq = max(mag(epsilonEq - epsilonEqOld)/mesh_.time().deltaT0().value(), 1e-3); 

                scalar term1 = KcoeffP[faceID] / pow(EiP[faceID], nP[faceID]) * pow(ddtEpsilonEq/1e-3, mP[faceID]);
                // scalar term1 = KcoeffP[faceID] / pow(EiP[faceID], nP[faceID]) * pow(1, mP[faceID]);
                scalar term2 = 1/(1-nP[faceID]);
                scalar sigmaYOld = KcoeffP[faceID]*pow(epsilonEqOld, nP[faceID]) * pow(ddtEpsilonEq/1e-3, mP[faceID]); 
                
                sigmaYP[faceID] = max(pow(term1, term2), sigmaYOld);
                HpP[faceID] = max(1/((1/nP[faceID]/KcoeffP[faceID])*pow(sigmaEq/KcoeffP[faceID],1/nP[faceID] - 1)*pow(1e-3/ddtEpsilonEq, mP[faceID]/nP[faceID]) - 1/EiP[faceID]), 0.0);
            }
        }
    } 
}

// ************************************************************************* //