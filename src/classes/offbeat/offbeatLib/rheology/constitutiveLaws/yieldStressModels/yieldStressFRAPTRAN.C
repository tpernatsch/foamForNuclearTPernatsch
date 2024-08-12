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

#include "yieldStressFRAPTRAN.H"
#include "addToRunTimeSelectionTable.H"

#include "fvm.H"
#include "fvc.H"
#include "PstreamReduceOps.H"
#include "Tuple2.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(yieldStressFRAPTRAN, 0);
    addToRunTimeSelectionTable
    (
        yieldStressModel, 
        yieldStressFRAPTRAN, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::yieldStressFRAPTRAN::yieldStressFRAPTRAN
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    yieldStressModel(mesh, lawDict),
    phiValue_(readScalar(lawDict.lookup("phi"))),
    CWValue_(readScalar(lawDict.lookup("CW"))),
    Kcoeff_
    (
        IOobject
        (
            "Kcoeff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("Kcoeff", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    m_
    (
        IOobject
        (
            "m",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("m", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    ),
    n_
    (
        IOobject
        (
            "n",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("n", dimless, 0),
        zeroGradientFvPatchField<scalar>::typeName
    )/*,
    EpsDot_
    (
        IOobject
        (
            "ddt(epsilon)",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fvc::ddt(mesh.lookupObject<symmTensorField>("epsilon"))
    )*/
{
    
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yieldStressFRAPTRAN::~yieldStressFRAPTRAN()
{}


// * * * * * * * * * * * * * * * * Selectors  * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::yieldStressFRAPTRAN::correctYieldStress
(
    const labelList& addr
)
{
    ///- Elastic modulus and temperature
    const volScalarField& E(mesh_.lookupObject<volScalarField>("E"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    //-Strain --> strain rate 
 //       Eps_(mesh.lookupObject<volSymmTensorField>"epsilon"),
  //  const symmTensorField& EpsDot(mesh_.lookupObject<volSymmTensorField>("ddt(epsilon)"));

/*    volScalarField EpsDot_XX(EpsDot.component(symmTensor::XX));
    volScalarField EpsDot_YY(EpsDot.component(symmTensor::YY));
    volScalarField EpsDot_ZZ(EpsDot.component(symmTensor::ZZ));*/


    phiValue_ *= 1e25;
    //- strain rate exponent m
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        const scalar T = Temp.internalField()[cellI];

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
    }

    //- strain hardening exponent n
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        const scalar T = Temp.internalField()[cellI];

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
    }

    //- Strength coefficient
    const scalar KCW = 0.546 * CWValue_;
    scalar Kphi = 0;
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        const scalar T = Temp.internalField()[cellI];

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

    }

    //- yieldStress update
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];
        const scalar Ei = E.internalField()[cellI];
 /*       const scalar EpsDotXXi = EpsDot_XX.internalField()[cellI];
        const scalar EpsDotYYi = EpsDot_YY.internalField()[cellI];
        const scalar EpsDotZZi = EpsDot_ZZ.internalField()[cellI];*/
        //- ddt(epsilon) is a symmTensor (obv). We choose the normal stress component which is the higher to maximize strain-rate effect...
        //- Doesnt change much in practice as the strain rate effect is small enough so that the difference on the yield stress in the XX, YY or ZZ direction is negligeable even if the strain rates in each directions are not stricly equals

//      const scalar EpsDotMax = std::max(EpsDotXXi, std::max(EpsDotYYi, EpsDotZZi));
        //scalar term1 = Kcoeff_[cellI] / pow(Ei, n_[cellI]) * pow(EpsDotMax, abs(m_[cellI]));
        scalar term1 = Kcoeff_[cellI] / pow(Ei, n_[cellI]) * pow(1, m_[cellI]);
        scalar term2 = 1/(1-n_[cellI]);
        
        sigmaY_[cellI] = pow(term1, term2);

        n_.correctBoundaryConditions();
        m_.correctBoundaryConditions();
        Kcoeff_.correctBoundaryConditions();

        const cell& c = mesh_.cells()[cellI];  

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

            if (patchID > -1 and sigmaY_.boundaryField()[patchID].size())
            {
                scalarField EiP = E.boundaryField()[patchID];
                scalarField TempP = Temp.boundaryField()[patchID];
                scalarField KcoeffP = Kcoeff_.boundaryField()[patchID];
                scalarField nP = n_.boundaryField()[patchID];
                scalarField mP = m_.boundaryField()[patchID];

                scalarField term1 = KcoeffP / pow(EiP, nP) * pow(1, mP);
                scalarField term2 = 1/(1-nP);
                scalarField sigmaY_valueP = pow(term1, term2);

                scalarField& sigmaYP(sigmaY_.boundaryFieldRef()[patchID]);

                const label faceID = mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                
                sigmaYP[faceID] = sigmaY_valueP[faceID];
            }
        }
    } 
}
// ************************************************************************* //