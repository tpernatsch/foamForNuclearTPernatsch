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

#include "ZircaloyLimback.H"
#include "addToRunTimeSelectionTable.H"
#include "globalFieldLists.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ZircaloyLimback, 0);
    addToRunTimeSelectionTable
    (
        creepModel, 
        ZircaloyLimback, 
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ZircaloyLimback::ZircaloyLimback
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    fastFluenceName_(lawDict.lookupOrDefault<word>("fastFluenceName", "fastFluence" )),
    fastFluxName_(lawDict.lookupOrDefault<word>("fastFluxName", "fastFlux" )),
    fastFluence_(nullptr),
    fastFlux_(nullptr),
    epsilonCreepIrr_(createOrLookup<symmTensor>(mesh, "epsilonCreepIrr")),
    epsilonCreepTh_(createOrLookup<symmTensor>(mesh, "epsilonCreepTh")),
    epsilonCreepPrim_(createOrLookup<symmTensor>(mesh, "epsilonCreepPrim")),
    epsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "epsilonCreepIrrEq")),
    DepsilonCreepIrrEq_(createOrLookup<scalar>(mesh, "DepsilonCreepIrrEq")),
    epsilonCreepThEq_(createOrLookup<scalar>(mesh, "epsilonCreepThEq")),
    DepsilonCreepThEq_(createOrLookup<scalar>(mesh, "DepsilonCreepThEq")),
    epsilonCreepPrimEq_(createOrLookup<scalar>(mesh, "epsilonCreepPrimEq")),
    DepsilonCreepPrimEq_(createOrLookup<scalar>(mesh, "DepsilonCreepPrimEq")),
    timeIndex_(0),
    relax_(lawDict.lookupOrDefault<scalar>("relax", 1.0 )),
    cladType_("SRA"),
    R(8.314),
    C0(3.557e-24),
    C1(0.85),
    C2(1.0),
    A(1.08e9),
    Q(201e3),
    n(2),
    A1(0.56),
    a(650),
    A3(1.3),
    A2(1.4e-27),
    C(52),
    b(0.109),
    d(-2.05),
    B(0.0216),
    D(35500)
{
    epsilonCreepIrr_.oldTime();
    epsilonCreepTh_.oldTime();
    epsilonCreepPrim_.oldTime();
    epsilonCreepIrrEq_.oldTime();
    DepsilonCreepIrrEq_.oldTime();
    epsilonCreepThEq_.oldTime();
    DepsilonCreepThEq_.oldTime();
    epsilonCreepPrimEq_.oldTime();
    DepsilonCreepPrimEq_.oldTime();

    if(lawDict.found("cladType"))
    {
        cladType_ = word(lawDict.lookup("cladType"));

        if (cladType_ == "SRA")
        {
            A = 1.08e9;
            Q = 201e3;
            n = 2.0;

            C0 = 3.557e-24;
            C1 = 0.85;
            C2 = 1;
        }
        else if (cladType_ == "RXA")
        {
            A = 5.47e8;
            Q = 198e3;
            n = 3.5;

            C0 = 1.654e-24;
            C1 = 0.85;
            C2 = 1;
        }
        else if (cladType_ == "PRA")
        {
            A = 7.06e8;
            Q = 199e3;
            n = 2.3;

            C0 = 2.714e-24;
            C1 = 0.85;
            C2 = 1;
        }
        else if (cladType_ == "ZIRLO")
        {
            A = 8.64e8;
            Q = 201e3;

            // Placeholder, for Zirlo n is a function of the stress
            n = 2.0;

            C0 = 2.846e-24;
            C1 = 0.85;
            C2 = 1;
        }
        else
        {
            FatalErrorInFunction
                << nl
                << "    Wrong \"cladType\" " << cladType_ << " for ZircaloyLimback."
                << "Available choices are:" << nl
                << "    - SRA i.e. stress relief annealed (Zr2 or Zr4), default choice" << nl
                << "    - RXA i.e. recrystallization annealed Zr2 or M5)" << nl
                << "    - PRA i.e. partially recrystallization annealed (Zr2)" << nl
                << "    - ZIRLO i.e. stress relief annealed ZIRLO" << nl
                << abort(FatalError);
        }
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ZircaloyLimback::~ZircaloyLimback()
{
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::ZircaloyLimback::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{    
    if(fastFluence_ == nullptr)
    {        
        if(mesh_.foundObject<volScalarField>(fastFluenceName_))
        {
            fastFluence_ = &mesh_.lookupObject<volScalarField>(fastFluenceName_);
        }
        else
        {
            FatalErrorInFunction
                << nl
                << "    Fast fluence field \"" << fastFluenceName_ << "\" requested by " 
                << "ZircaloyLimback but not found in registry." << nl
                << abort(FatalError);
        }

        if(mesh_.foundObject<volScalarField>(fastFluxName_))
        {
            fastFlux_ = &mesh_.lookupObject<volScalarField>(fastFluxName_);
        }
        else
        {
            FatalErrorInFunction
                << nl
                << "    Fast flux field \"" << fastFluxName_ << "\" requested by " 
                << "ZircaloyLimback but not found in registry." << nl
                << abort(FatalError);
        }
    }   

    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    scalar tol(1.0);

    //- Irradiation time step [hr]
    scalar deltaT(0);
    
    if(mesh_.time().timeIndex() > 0)
    {
        deltaT = mesh_.time().deltaTValue()/3600.0;

        const scalarField& epsilonCreepIrrEqOld_ = 
        epsilonCreepIrrEq_.oldTime().internalField();
        const scalarField& epsilonCreepThEqOld_ = 
        epsilonCreepThEq_.oldTime().internalField();
        const scalarField& epsilonCreepPrimEqOld_ = 
        epsilonCreepPrimEq_.oldTime().internalField();

        const scalarField& fastFluenceI_ = fastFluence_->internalField();
        const scalarField& fastFluxI_ = fastFlux_->internalField();

        const scalarField& epsilonCreepEqOld_ = 
        epsilonCreepEq_.oldTime().internalField();

        const symmTensorField& epsilonCreepIrrOldI = 
        epsilonCreepIrr_.oldTime().internalField();
        const symmTensorField& epsilonCreepThOldI = 
        epsilonCreepTh_.oldTime().internalField();
        const symmTensorField& epsilonCreepPrimOldI = 
        epsilonCreepPrim_.oldTime().internalField();
        const symmTensorField& epsilonCreepOldI = 
        epsilonCreep_.oldTime().internalField();
    
        //- Adjust elastic strain tensor by subtracting old creep values    
        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];

            epsilonEl[cellI] -= epsilonCreepOldI[cellI];

            const cell& c = mesh_.cells()[cellI];     

            forAll(c, faceI)
            {
                const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                {
                    if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                    {
                        const label faceID = 
                        mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                        
                        symmTensorField& epsilonElP = 
                        epsilonEl.boundaryFieldRef()[patchID];
                        const symmTensorField& epsilonCreepOldP = 
                        epsilonCreep_.oldTime().boundaryField()[patchID];

                        epsilonElP[faceID] -= epsilonCreepOldP[faceID];
                    }
                }   
            }
        }    

        if(mesh_.time().timeIndex() > timeIndex_)
        { 
            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];
                DepsilonCreepEq_.ref()[cellI] *= 0;
                DepsilonCreepIrrEq_.ref()[cellI] *= 0;
                DepsilonCreepThEq_.ref()[cellI] *= 0;
                DepsilonCreepPrimEq_.ref()[cellI] *= 0;

                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                    {
                        if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                            scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepIrrEqP = DepsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepThP = DepsilonCreepThEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepPrimEqP = DepsilonCreepPrimEq_.boundaryFieldRef()[patchID];

                            DepsilonCreepEqP[faceID] *= 0;
                            DepsilonCreepIrrEqP[faceID] *= 0;
                            DepsilonCreepThP[faceID] *= 0;
                            DepsilonCreepPrimEqP[faceID] *= 0;
                        }
                    }
                }
            }

            // DepsilonCreep_() *= 0;

            timeIndex_ = mesh_.time().timeIndex();
        }

        do
        {
            scalar maxRate(0);
            scalar maxFlux(0);
            scalar maxTemp(0);
            scalar maxStress(0);
            // scalar volTot(0);

            forAll(addr, addrI)
            {
                const label cellI = addr[addrI];

                const symmTensor sigmaI = 
                2*mu[cellI]*epsilonEl[cellI] + lambda[cellI]*tr(epsilonEl[cellI])*I;

                //- Calculate the von Mises stress [Pa]
                const scalar sigmaEff =
                min
                (
                    sqrt((3.0/2.0)
                    *magSqr(dev(sigmaI))) 
                    - 3*mu[cellI]*DepsilonCreepEq_.internalField()[cellI] 
                    , 3e10
                );
                
                maxStress = max(sigmaEff, maxStress);

                const scalar DepsilonCreepEqPrev = DepsilonCreepEq_.internalField()[cellI];

                //- Reference to temperature field [K]
                const scalar T = min(Temp.internalField()[cellI], 700.0);

                maxTemp = max(T, maxTemp);

            /****** Calculate irradiation creep rate [1/s]  *********/  

                //- Fast flux 
                scalar phi = fastFluxI_[cellI]; 
                maxFlux = max(phi, maxFlux);

                const scalar creepIrrRate = C0*pow(phi*1e4, C1)*pow(sigmaEff/1e6, C2);

                DepsilonCreepIrrEq_.ref()[cellI] =  
                (1 - relax_)*DepsilonCreepIrrEq_.internalField()[cellI] + relax_*creepIrrRate*deltaT;

            /****** Calculate secondary thermal creep rate [1/hr] ********/

                const scalar E = 1.148e5 - 59.9*T;

                const scalar ai = a*(1 - A1*(1 - exp(-A2*pow(fastFluenceI_[cellI],A3))));

                scalar creepThRate(0);

                // Correct n for zirlo claddings
                if(cladType_ == "ZIRLO")
                {
                    if(sigmaEff/1e6 < 220)
                    {
                        n = 2.0;
                    }
                    else if(sigmaEff/1e6 < 400)
                    {
                        n = 2.6;
                    }
                    else
                    {
                        n = 1.2667 + 3.333e-3*sigmaEff/1e6;
                    }
                }

                creepThRate = A*E/T*pow(Foam::sinh(ai*sigmaEff/1e6/E), n)*exp(-Q/R/T);

                DepsilonCreepThEq_.ref()[cellI] =   
                (1 - relax_)*DepsilonCreepThEq_.internalField()[cellI] + relax_*creepThRate*deltaT;

            /****** Calculate total secondary creep rate *********/
                const scalar creepSecRate = creepIrrRate + creepThRate;

            /****** Calculate primary creep *********/

                if(creepSecRate > VSMALL)
                {
                    //- Saturated primary creep
                    const scalar epsilonPrim_sat = 
                    B*pow(creepSecRate, b)*pow(2 - Foam::tanh(D*creepSecRate), d);

                    if(epsilonCreepPrimEqOld_[cellI] < (epsilonPrim_sat-SMALL))
                    {
                        //- Time constant
                        const scalar tau = 
                        1/C/C*pow
                        (
                            Foam::log
                            (
                                1 - epsilonCreepPrimEqOld_[cellI]/epsilonPrim_sat
                            )
                            ,2
                        )
                        *1/creepSecRate;

                        //- Primary creep
                        DepsilonCreepPrimEq_.ref()[cellI] = 
                        (1 - relax_)*DepsilonCreepPrimEq_.internalField()[cellI]
                        + relax_*
                        (
                            epsilonPrim_sat*(1 - exp(-C*sqrt(creepSecRate*(tau + deltaT)))) 
                            - epsilonCreepPrimEqOld_[cellI]
                        );
                        
                        epsilonCreepPrimEq_.ref()[cellI] = epsilonCreepPrimEqOld_[cellI] + DepsilonCreepPrimEq_.internalField()[cellI];
                    }
                }

            /******* Total creep ********/
                // const scalar alpha = 
                // max(1 - epsilonCreepPrimEq_.ref()[cellI]/epsilonPrim_sat, SMALL);

                // const scalar alphaOld = 
                // 1 - epsilonCreepPrimEqOld_[cellI]/epsilonPrim_sat;

                const scalar fTrial = 
                creepSecRate + DepsilonCreepPrimEq_.internalField()[cellI]/deltaT;

                // const scalar fTrialDev = (C2*creepIrrRate/sigmaEff*1e6 + n*ai*creepThRate/E*(1.0/Foam::tanh(ai*sigmaEff/1e6/E)) )*
                //                          (1 + 1/deltaT*(B/creepSecRate - (d*D*pow(Foam::cosh(D*creepSecRate), -2.0))/(2 - Foam::tanh(D*creepSecRate)))*
                //                          (epsilonCreepPrimEq_.ref()[cellI] - alpha*Foam::log(alphaOld)*epsilonCreepPrimEqOld_[cellI]/(alphaOld*Foam::log(alpha))) - epsilonPrim_sat*C*C*alpha/Foam::log(alpha));

                //DepsilonCreepEq_.ref()[cellI] = (1 - relax_)*DepsilonCreepEq_.ref()[cellI] + relax_*(DepsilonCreepEq_.ref()[cellI] - (DepsilonCreepEq_.ref()[cellI] - deltaT*fTrial)/(1 + 3*mu[cellI]/1e6*fTrialDev*deltaT));
                
                DepsilonCreepEq_.ref()[cellI] = 
                (1 - relax_)*DepsilonCreepEq_.internalField()[cellI] + relax_*fTrial*deltaT;
                
                maxRate = max(fTrial, maxRate);

                epsilonCreepIrrEq_.ref()[cellI] = 
                epsilonCreepIrrEqOld_[cellI] + DepsilonCreepIrrEq_.internalField()[cellI];

                epsilonCreepIrr_.ref()[cellI] = 
                epsilonCreepIrrOldI[cellI] + 1.5*DepsilonCreepIrrEq_.internalField()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.internalField()[cellI])/max(sigmaEff, SMALL);


                epsilonCreepThEq_.ref()[cellI] = 
                epsilonCreepThEqOld_[cellI] + DepsilonCreepThEq_.internalField()[cellI];

                epsilonCreepTh_.ref()[cellI] = 
                epsilonCreepThOldI[cellI] + 1.5*DepsilonCreepThEq_.internalField()[cellI]*dev(sigmaI  
                    - 2*mu[cellI]*DepsilonCreep_.internalField()[cellI])/max(sigmaEff, SMALL); 



                epsilonCreepPrim_.ref()[cellI] = 
                epsilonCreepPrimOldI[cellI] + 1.5*DepsilonCreepPrimEq_.internalField()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.internalField()[cellI])/max(sigmaEff, SMALL);


                epsilonCreepEq_.ref()[cellI] = 
                epsilonCreepEqOld_[cellI] + DepsilonCreepEq_.internalField()[cellI];

                DepsilonCreep_.ref()[cellI] = 
                (1 - relax_)*DepsilonCreep_.internalField()[cellI] + relax_*1.5*DepsilonCreepEq_.internalField()[cellI]*dev(sigmaI 
                    - 2*mu[cellI]*DepsilonCreep_.internalField()[cellI])/max(sigmaEff, SMALL); 

                epsilonCreep_.ref()[cellI] = 
                epsilonCreepOldI[cellI] + DepsilonCreep_.internalField()[cellI]; 


                const scalar tolTr = min(1, mag( (DepsilonCreepEq_.internalField()[cellI] - DepsilonCreepEqPrev)
                    / (max(DepsilonCreepEqPrev , SMALL))));

                tol = min(tol, tolTr );

                const cell& c = mesh_.cells()[cellI];  

                forAll(c, faceI)
                {
                    const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);

                    {
                        if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                        {
                            const label faceID = 
                            mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);

                            const scalarField& Tp = Temp.boundaryField()[patchID];
                            const scalarField& muP = mu.boundaryField()[patchID];
                            const scalarField& lambdaP = lambda.boundaryField()[patchID];
                            const scalarField& fastFluenceP = fastFluence_->boundaryField()[patchID];
                            const scalarField& fastFluxP = fastFlux_->boundaryField()[patchID];
                            const symmTensorField& epsilonElP = epsilonEl.boundaryField()[patchID];

                            symmTensorField& epsilonCreepP = epsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& DepsilonCreepP = DepsilonCreep_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepIrrP = epsilonCreepIrr_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepThP = epsilonCreepTh_.boundaryFieldRef()[patchID];
                            symmTensorField& epsilonCreepPrimP = epsilonCreepPrim_.boundaryFieldRef()[patchID];

                            scalarField& epsilonCreepEqP = epsilonCreepEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepIrrEqP = epsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepThEqP = epsilonCreepThEq_.boundaryFieldRef()[patchID];
                            scalarField& epsilonCreepPrimEqP = epsilonCreepPrimEq_.boundaryFieldRef()[patchID];

                            const scalarField& epsilonCreepIrrEqOldP = epsilonCreepIrrEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepThEqOldP = epsilonCreepThEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepPrimEqOldP = epsilonCreepPrimEq_.oldTime().boundaryField()[patchID];
                            const scalarField& epsilonCreepEqOldP = epsilonCreepEq_.oldTime().boundaryField()[patchID];

                            const symmTensorField& epsilonCreepIrrOldP = epsilonCreepIrr_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepThOldP = epsilonCreepTh_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepPrimOldP = epsilonCreepPrim_.oldTime().boundaryField()[patchID];
                            const symmTensorField& epsilonCreepOldP = epsilonCreep_.oldTime().boundaryField()[patchID];

                            scalarField& DepsilonCreepIrrEqP = DepsilonCreepIrrEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepThEqP = DepsilonCreepThEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepPrimEqP = DepsilonCreepPrimEq_.boundaryFieldRef()[patchID];
                            scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];

                            //- Calculate the von Mises stress [Pa]
                            const symmTensor sigmaPI = 
                            2*muP[faceID]*epsilonElP[faceID] + lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                            const scalar sigmaEff = 
                            min(sqrt((3.0/2.0)*magSqr(dev(sigmaPI))) - 3*muP[faceID]*DepsilonCreepEqP[faceID], 3e10);

                            const scalar DepsilonCreepEqPPrev = DepsilonCreepEqP[faceID];

                            //- Reference to temperature field [K]
                            const scalar T = min(Tp[faceID], 700.0);

                        /****** Calculate irradiation creep rate [1/s]  *********/   

                            //- Fast flux 
                            scalar phi = fastFluxP[faceID];

                            const scalar creepIrrRate = C0*pow(phi*1e4, C1)*pow(sigmaEff/1e6, C2);

                            DepsilonCreepIrrEqP[faceID] =  (1 - relax_)*DepsilonCreepIrrEqP[faceID] + relax_*creepIrrRate*deltaT;

                        /****** Calculate secondary thermal creep rate [1/hr] ********/

                            const scalar E = 1.148e5 - 59.9*T;

                            const scalar ai = a*(1 - A1*(1 - exp(-A2*pow(fastFluenceP[faceID],A3))));
                            //const scalar ai = 650;

                            scalar creepThRate(0);

                            // Correct n for zirlo claddings
                            if(cladType_ == "ZIRLO")
                            {
                                if(sigmaEff/1e6 < 220)
                                {
                                    n = 2.0;
                                }
                                else if(sigmaEff/1e6 < 400)
                                {
                                    n = 2.6;
                                }
                                else
                                {
                                    n = 1.2667 + 3.333e-3*sigmaEff/1e6;
                                }
                            }

                            creepThRate = A*E/T*pow(Foam::sinh(ai*sigmaEff/1e6/E), n)*exp(-Q/R/T);

                            DepsilonCreepThEqP[faceID] = (1 - relax_)*DepsilonCreepThEqP[faceID] + relax_*creepThRate*deltaT;

                        /****** Calculate total secondary creep rate *********/

                            const scalar creepSecRate = creepIrrRate + creepThRate;

                        /****** Calculate primary creep *********/

                            if(creepSecRate > VSMALL)
                            {

                                //- Saturated primary creep
                                const scalar epsilonPrim_sat = B*pow(creepSecRate, b)*pow(2 - Foam::tanh(D*creepSecRate), d);

                                if(epsilonCreepPrimEqOldP[faceID] < (epsilonPrim_sat-SMALL))
                                {
                                    //- Time constant
                                    const scalar tau = 1/C/C*pow(Foam::log(1 - epsilonCreepPrimEqOldP[faceID]/epsilonPrim_sat), 2)*1/creepSecRate;

                                    //- Primary creep
                                    DepsilonCreepPrimEqP[faceID] = 
                                    (1 - relax_)*DepsilonCreepPrimEqP[faceID] + relax_*(epsilonPrim_sat*(1 - exp(-C*sqrt(creepSecRate*(tau + deltaT)))) - epsilonCreepPrimEqOldP[faceID]);
                                

                                    epsilonCreepPrimEqP[faceID] = 
                                    DepsilonCreepPrimEqP[faceID] + epsilonCreepPrimEqOldP[faceID] ;                        
                                }
                            }

                        /******* Total creep ********/
                                    // const scalar alpha = 
                                    // max(1 - epsilonCreepPrimEqP[faceID]/epsilonPrim_sat, SMALL);

                                    // const scalar alphaOld = 
                                    // 1 - epsilonCreepPrimEqOldP[faceID]/epsilonPrim_sat;

                            const scalar fTrial = 
                            creepSecRate + DepsilonCreepPrimEqP[faceID]/deltaT;

                            // const scalar fTrialDev = (C2*creepIrrRate/sigmaEff*1e6 + n*ai*creepThRate/E*(1.0/Foam::tanh(ai*sigmaEff/1e6/E)) )*
                            //                          (1 + 1/deltaT*(B/creepSecRate - (d*D*pow(Foam::cosh(D*creepSecRate), -2.0))/(2 - Foam::tanh(D*creepSecRate)))*
                            //                          (epsilonCreepPrimEqP[faceID] - alpha*Foam::log(alphaOld)*epsilonCreepPrimEqOldP[faceID]/(alphaOld*Foam::log(alpha))) - epsilonPrim_sat*C*C*alpha/Foam::log(alpha));

                            //DepsilonCreepEqP[faceID] = (1 - relax_)*DepsilonCreepEqP[faceID] + relax_*(DepsilonCreepEqP[faceID] - (DepsilonCreepEqP[faceID] - deltaT*fTrial)/(1 + 3*muP[faceID]/1e6*fTrialDev*deltaT));
                            DepsilonCreepEqP[faceID] = (1-relax_)*DepsilonCreepEqP[faceID] + relax_*fTrial*deltaT;


                            epsilonCreepIrrEqP[faceID] = 
                            DepsilonCreepIrrEqP[faceID] + epsilonCreepIrrEqOldP[faceID] ;
                            
                            epsilonCreepIrrP[faceID] = 
                            epsilonCreepIrrOldP[faceID] + 1.5*DepsilonCreepIrrEqP[faceID]*dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);


                            epsilonCreepThEqP[faceID] = 
                            DepsilonCreepThEqP[faceID] + epsilonCreepThEqOldP[faceID] ;
                            
                            epsilonCreepThP[faceID] = 
                            epsilonCreepThOldP[faceID] + 1.5*DepsilonCreepThEqP[faceID]*dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);

                            
                            epsilonCreepPrimP[faceID] = 
                            epsilonCreepPrimOldP[faceID] + 1.5*DepsilonCreepPrimEqP[faceID]*dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL);


                            epsilonCreepEqP[faceID] = 
                            DepsilonCreepEqP[faceID] + epsilonCreepEqOldP[faceID] ;

                            DepsilonCreepP[faceID] = 
                            (1 - relax_)*DepsilonCreepP[faceID] + relax_*1.5*DepsilonCreepEqP[faceID]*dev(sigmaPI - 2*muP[faceID]*DepsilonCreepP[faceID])/max(sigmaEff, SMALL) ;
                            
                            epsilonCreepP[faceID] = 
                            epsilonCreepOldP[faceID] + DepsilonCreepP[faceID];

                            
                            const scalar tolTr = min(1, mag( (DepsilonCreepEqP[faceID] - DepsilonCreepEqPPrev)/ (max(DepsilonCreepEqPPrev , SMALL))));

                            tol = min(tol, tolTr );
                        }
                    }
                }
            }

        }while(tol > 1);
    }    
    
    forAll(addr, addrI)
    {
        const label cellI = addr[addrI];

        epsilonEl[cellI] -= DepsilonCreep_.internalField()[cellI];     
        
        const cell& c = mesh_.cells()[cellI];         

        forAll(c, faceI)
        {
            const label patchID = mesh_.boundaryMesh().whichPatch(c[faceI]);   

            {
                if (patchID > -1 and epsilonEl.boundaryField()[patchID].size())
                {
                    const label faceID = 
                    mesh_.boundaryMesh()[patchID].whichFace(c[faceI]);
                    
                    symmTensorField& epsilonElP = 
                    epsilonEl.boundaryFieldRef()[patchID];
                    const symmTensorField& DepsilonCreepP = 
                    DepsilonCreep_.boundaryField()[patchID];

                    epsilonElP[faceID] -= DepsilonCreepP[faceID];
                }
            }   
        }
    }    
}

// ************************************************************************* //


