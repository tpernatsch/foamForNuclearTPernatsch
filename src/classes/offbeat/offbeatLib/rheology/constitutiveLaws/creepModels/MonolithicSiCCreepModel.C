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

#include "MonolithicSiCCreepModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MonolithicSiCCreepModel, 0);
    addToRunTimeSelectionTable
    (
        creepModel,
        MonolithicSiCCreepModel,
        dictionary
    );
}

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MonolithicSiCCreepModel::MonolithicSiCCreepModel
(
    const fvMesh& mesh,
    const dictionary& lawDict
)
:
    creepModel(mesh, lawDict),
    timeIndex_(0),
    relax_(lawDict.lookupOrDefault<scalar>("relaxCreep", 1.0 )),
    fastFluxName_(lawDict.lookupOrDefault<word>("fastFluxName", "fastFlux" )),
    K_(lawDict.lookupOrDefault<scalar>("creepCoefficient", -1)),
    vc_(lawDict.lookupOrDefault<scalar>("poissonRationForCreep", 0.5)),
    fcf_(lawDict.lookupOrDefault<scalar>("fluxConversionFactor", 1.0)),
    par1(2.38),
    par2(1.9),
    par3(2.193e-29),
    par4(4.85e-32),
    par5(4.0147e-35)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::MonolithicSiCCreepModel::~MonolithicSiCCreepModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::MonolithicSiCCreepModel::correctCreep
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

            //- Calculate the von Mises stress [Pa]
            const scalar sigmaEff =
            min
            (
                sqrt((3.0/2.0)
                *magSqr(dev(sigmaI)))
                - 3*mu[cellI]*DepsilonCreepEq_.internalField()[cellI]
                , 3e10
            );

            const scalar DepsilonCreepEqPrev = DepsilonCreepEq_.ref()[cellI];

            //- Reference to temperature field [Celsius]. The remperature here only
            //- affect creep coefficient. It is bounded at 640.0 C and 900.0 C.
            const scalar T = max(min(Temp.internalField()[cellI]-273.15, 900.0),640.0);

            maxTemp = max(T, maxTemp);

            // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #
            scalar K;

            if (K_ < 0.0) //The case where the creep coefficient are not provided by the user
            {
              K = (4e-37-2e-37)/(900.0-640.0)*(T-640.0)+2e-37;
            }
            else
            {
              K = K_;
            }

            //- Calculate creep strain rate
            const scalar creepRate = K*sigmaEff*phi[cellI]*1e4*fcf_;

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

                        //- Calculate the von Mises stress [Pa]
                        const scalar sigmaEff =
                        min
                        (
                            sqrt((3.0/2.0)
                            *magSqr(dev(sigmaI)))
                            - 3*mu[cellI]*DepsilonCreepEq_.internalField()[cellI]
                            , 3e10
                        );

                        const scalar DepsilonCreepEqPPrev = DepsilonCreepEqP[faceID];
                        // * # * # * # * # * # * # * # * # * # * # * # * # * # * # * # * #
                        //- Reference to temperature field [Celsius]
                        const scalar T = max(min(Tp[faceID]-273.15, 900.0), 640.0);

                        scalar Kp;

                        if (K_ < 0.0) //The case where the creep coefficient are not provided by the user
                        {
                          Kp = (4e-37-2e-37)/(900.0-640.0)*(T-640.0)+2e-37;
                        }
                        else
                        {
                          Kp = K_;
                        }

                        //- Calculate creep strain rate with power law
                        const scalar creepRate = Kp*sigmaEff*phiP[faceID]*1e4*fcf_;

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

// ************************************************************************* //
