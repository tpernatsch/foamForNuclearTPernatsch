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

#include "PARFUMEBufferCreepModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PARFUMEBufferCreepModel, 0);
    addToRunTimeSelectionTable
    (
        creepModel,
        PARFUMEBufferCreepModel,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PARFUMEBufferCreepModel::PARFUMEBufferCreepModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    timeIndex_(0),
    relax_(lawDict.lookupOrDefault<scalar>("relaxCreep", 1.0 )),
    densityName_(lawDict.lookupOrDefault<word>("densityName", "rho" )),
    fastFluxName_(lawDict.lookupOrDefault<word>("fastFluxName", "fastFlux" )),
    C_amp(lawDict.lookupOrDefault<scalar>("creepAmplificationCoefficient", 2.0)),
    vc_(lawDict.lookupOrDefault<scalar>("poissonRationForCreep", 0.5)),
    fcf_(lawDict.lookupOrDefault<scalar>("fluxConversionFactor", 1.0)),
    par1(2.38),
    par2(1.9),
    par3(2.193e-29),
    par4(4.85e-32),
    par5(4.0147e-35)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::PARFUMEBufferCreepModel::~PARFUMEBufferCreepModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::PARFUMEBufferCreepModel::correctCreep
(
    volSymmTensorField& epsilonEl,
    const labelList& addr
)
{
    const volScalarField& mu(mesh_.lookupObject<volScalarField>("mu"));
    const volScalarField& lambda(mesh_.lookupObject<volScalarField>("lambda"));
    const volScalarField& Temp(mesh_.lookupObject<volScalarField>("T"));

    scalar tol(1.0);
    //- Irradiation time step [s]
    scalar deltaT(0);

    if(mesh_.time().timeIndex() > 0)
    {
        const volScalarField& rho = mesh_.lookupObject<volScalarField>(densityName_);
        const volScalarField& phi = mesh_.lookupObject<volScalarField>(fastFluxName_);

        deltaT = mesh_.time().deltaTValue();

        const scalarField& epsilonCreepEqOld_ = epsilonCreepEq_.oldTime().internalField();
        const symmTensorField& epsilonCreepOldI = epsilonCreep_.oldTime().internalField();

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

                            DepsilonCreepEqP[faceID] *= 0;
                        }
                    }
                }
            }

            timeIndex_ = mesh_.time().timeIndex();
        }

        // Do-while loop to correct creep strain
        do
        {
        scalar maxTemp(0);
        //scalar maxStress(0);
        //scalar volTot(0);

        forAll(addr, addrI)
        {
            const label cellI = addr[addrI];

            const symmTensor sigmaI = 2*mu[cellI]*epsilonEl[cellI] + lambda[cellI]*tr(epsilonEl[cellI])*I;

            //- Calculate the principal stress [Pa]
            scalar p1, p2, p3;
            bool s;
            computePrincipalStress(sigmaI, p1, p2, p3, s);
            //- Calculate the effective stress;
            //- Formula taken from "D. Petti, P. Martin, M. Phelip, and R. Ballinger.
            //- Development of improved models and designs for coated-particle gas
            //- reactor fuels. Technical Report INL/EXT-05-02615, Idaho National
            //- Laboratory, December 2004."
            const scalar sigmaEff = min(p1-vc_*(p2+p3)-3*mu[cellI]*DepsilonCreepEq_.ref()[cellI], 3e10);

            const scalar DepsilonCreepEqPrev = DepsilonCreepEq_.ref()[cellI];

            //- Reference to temperature field [Celsius]
            const scalar T = Temp.internalField()[cellI]-273.15;

            maxTemp = max(T, maxTemp);

            // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #
            //- Irradiation-induced creep coefficient
            const  scalar Ks = C_amp*(1+par1*(par2-rho[cellI]/1e3))*(par3-par4*T+par5*T*T);
            //- Calculate creep strain rate
            const scalar creepRate = Ks*sigmaEff*1e-6*phi[cellI]*1e4*fcf_;

            //- Calculate creep strain
            DepsilonCreepEq_.ref()[cellI] =
            (1-relax_)*DepsilonCreepEq_.ref()[cellI] + relax_*creepRate*deltaT;

            // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

            epsilonCreepEq_.ref()[cellI] =
            epsilonCreepEqOld_[cellI] + DepsilonCreepEq_.ref()[cellI];

            DepsilonCreep_.ref()[cellI] =
            (1 - relax_)*DepsilonCreep_.ref()[cellI] + relax_*1.5*DepsilonCreepEq_.ref()[cellI]*dev(sigmaI
                - 2*mu[cellI]*DepsilonCreep_.ref()[cellI])/max(sigmaEff, SMALL);

            epsilonCreep_.ref()[cellI] =
            epsilonCreepOldI[cellI] + DepsilonCreep_.ref()[cellI];


            const scalar tolTr = min(1, mag( (DepsilonCreepEq_.ref()[cellI] - DepsilonCreepEqPrev)
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
                        const scalarField& rhoP = rho.boundaryField()[patchID];
                        const scalarField& phiP = phi.boundaryField()[patchID];
                        const symmTensorField& epsilonElP = epsilonEl.boundaryField()[patchID];

                        symmTensorField& epsilonCreepP = epsilonCreep_.boundaryFieldRef()[patchID];
                        symmTensorField& DepsilonCreepP = DepsilonCreep_.boundaryFieldRef()[patchID];

                        scalarField& epsilonCreepEqP = epsilonCreepEq_.boundaryFieldRef()[patchID];

                        const scalarField& epsilonCreepEqOldP = epsilonCreepEq_.oldTime().boundaryField()[patchID];

                        const symmTensorField& epsilonCreepOldP = epsilonCreep_.oldTime().boundaryField()[patchID];

                        scalarField& DepsilonCreepEqP = DepsilonCreepEq_.boundaryFieldRef()[patchID];

                        const symmTensor sigmaPI =
                        2*muP[faceID]*epsilonElP[faceID] + lambdaP[faceID]*tr(epsilonElP[faceID])*I;

                        //- Calculate the principal stress [Pa]
                        scalar p1p, p2p, p3p;
                        bool sp;
                        computePrincipalStress(sigmaPI, p1p, p2p, p3p, sp);

                        const scalar sigmaEff = min(p1p-vc_*(p2p+p3p)-3*muP[faceID]*DepsilonCreepEqP[faceID], 3e10);

                        const scalar DepsilonCreepEqPPrev = DepsilonCreepEqP[faceID];
                        // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #
                        //- Reference to temperature field [Celsius]
                        const scalar T = Tp[faceID] - 273.15;
                        //- Irradiation-induced creep coefficient
                        const scalar Ks = C_amp*(1+par1*(par2-rhoP[faceID]/1e3))*(par3-par4*T+par5*T*T);

                        //- Calculate creep strain rate with power law
                        const scalar creepRate = Ks*sigmaEff*1e-6*phiP[faceID]*1e4*fcf_;

                        //- Calculate creep strain
                        DepsilonCreepEqP[faceID] =
                        (1-relax_)*DepsilonCreepEqP[faceID] + relax_*creepRate*deltaT;

                        // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #

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

        epsilonEl[cellI] -= DepsilonCreep_[cellI];

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
    timeIndex_ = mesh_.time().timeIndex();
}

//TO DO: The computation of principal stress should be coded elsewhere so that
//it does not need to be written again and again in different source files.
void Foam::PARFUMEBufferCreepModel::computePrincipalStress
(
    symmTensor s, //stress tensor
    scalar& stressPrinciple1,
    scalar& stressPrinciple2,
    scalar& stressPrinciple3,
    bool& solved //whether the principle can be solved, for now it is not be used.
)
{
  //- Parameters for the cubic equation SP³-I1*SP²+I2*SP-I3=0 (SP: Stress Principal)
  scalar I1 = s.xx()+s.yy()+s.zz();
  scalar I2(s.xx()*s.yy()+s.xx()*s.zz()+s.yy()*s.zz()
            -s.xy()*s.xy()-s.yz()*s.yz()-s.xz()*s.xz());
  scalar I3(det(s));


  //--------------------------Solve the cubic equation--------------------------
  //- Convert the cubic equation into the form y³+py+q=0, where y = SP-I1/3;
  scalar p((3.0*I2-I1*I1)/3.0);
  scalar q((9.0*I1*I2-27.0*I3-2.0*pow(I1, 3))/27.0);

  scalar a = pow(q/2.0, 2);
  scalar b = pow(p/3.0, 3);
  scalar delta = a + b;

  bool noNeedtoSolve = false;
  scalar min_ab = max(min(mag(a), mag(b)), SMALL);

  const scalar minNstress = min(min(mag(s.xx()), mag(s.yy())), mag(s.zz()));
  const scalar maxSstress = max(max(mag(s.xy()), mag(s.xz())), mag(s.yz()));

  //- For the case when the stress tensor is already in the form of principal stress.
  //- Then there's no need to do the calculation, if not, there will cause errors
  //- because of the calculation error
  if (mag(delta/min_ab) < 1e-5 || maxSstress/minNstress < 1e-3)
  {
    noNeedtoSolve = true;
    delta = -1; //Just need to give a negative value

  }


  scalar SP1, SP2, SP3;


  //- delta must < 0 to have 3 different real roots; However, at the beginning of
  //- the iteration, the stress tensor might not be realistic (true?) so that the principal
  //- could not be able to be solved. In this case, we do not calculate the creep(stress = 0.0).
  //- Or, maybe we can use Von Mises stress instead before the principal stress
  //- can be solved? (TO DO)
  if (delta > 0.0)
  {
    stressPrinciple1 = 0.0;
    stressPrinciple2 = 0.0;
    stressPrinciple3 = 0.0;
    solved = false;

    return;
  }

  if (noNeedtoSolve == true)
  {
    SP1 = s.xx();
    SP2 = s.yy();
    SP3 = s.zz();
  }
  else
  {
    scalar r = max(sqrt(-pow(p/3.0,3)), SMALL);

    scalar theta = 1.0/3.0*acos(-q/2.0/r);
    //- Three real roots
    scalar y1(2*pow(r,1.0/3.0)*cos(theta));
    scalar y2(2*pow(r,1.0/3.0)*cos(theta+2.0/3.0*M_PI));
    scalar y3(2*pow(r,1.0/3.0)*cos(theta+4.0/3.0*M_PI));
    //- Three principle stress
    SP1 = y1+I1/3.0;
    SP2 = y2+I1/3.0;
    SP3 = y3+I1/3.0;
  }

  //------------ Find the max, mid, min stress; SP1 > SP2 > SP3----------------
  scalar t;
  if (SP3 > SP2)
  {
    t = SP3;
    SP3 = SP2;
    SP2 = t;
  }
  if (SP3 > SP1)
  {
    t = SP3;
    SP3 = SP1;
    SP1 = t;
  }
  if (SP2 > SP1)
  {
    t = SP2;
    SP2 = SP1;
    SP1 = t;
  }
  //- The first(max) principle stress
  stressPrinciple1 = SP1;
  //- The second(mid) principle stress
  stressPrinciple2 = SP2;
  //- The third(min) principle stress
  stressPrinciple3 = SP3;
  solved = true;
}

// ************************************************************************* //
