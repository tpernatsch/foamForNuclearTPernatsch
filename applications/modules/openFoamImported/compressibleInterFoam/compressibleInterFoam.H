/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2512                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2025 OpenCFD Ltd.          |
-------------------------------------------------------------------------------
License
    This file is part of GeN-Foam.

    GeN-Foam is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    GeN-Foam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    This offering is not approved or endorsed by the OpenFOAM Foundation nor
    OpenCFD Limited, producer and distributor of the OpenFOAM(R)software via
    www.openfoam.com, and owner of the OPENFOAM(R) and OpenCFD(R) trademarks.

    This particular snippet of code is developed according to the developer's
    knowledge and experience in OpenFOAM. The users should be aware that
    there is a chance of bugs in the code, though we've thoroughly test it.
    The source code may not be in the OpenFOAM coding style, and it might not
    be making use of inheritance of classes to full extent.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::thermalHydraulicsModels::compressibleInterFoam

Description
    This class encapsulates a single fluid + 1 stationary structure
    Eulerian solver with a user-defineable regime-map based phyiscs modelling.

    Author:
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    compressibleInterFoam.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleInterFoam_H
#define compressibleInterFoam_H

#include "fvCFD.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "compressibleInterPhaseTransportModel.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"
#include "solver.H"
#include "IOMRFZoneList.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
//class regime;
//class regimeMapModel;

namespace solvers
{

/*---------------------------------------------------------------------------*\
                      Class compressibleInterFoam Declaration
\*---------------------------------------------------------------------------*/

class compressibleInterFoam
:
    public solver
{
private:

    pimpleControl pimple_;

    //- Phases

    volScalarField p_rgh_;
    volVectorField U_;
    surfaceScalarField phi_;
    twoPhaseMixtureThermo mixture_;
    volScalarField& alpha1_;
    volScalarField& alpha2_;
    const volScalarField& rho1_;
    const volScalarField& rho2_;
    volScalarField rho_;
    dimensionedScalar pMin_;

    //- Gravity related fields
    uniformDimensionedVectorField g_;
    uniformDimensionedScalarField hRef_;
    volScalarField gh_;
    surfaceScalarField ghf_;

    surfaceScalarField rhoPhi_;
    volScalarField dgdt_;

    // //Create alphaFluxes
    // IOobject alphaPhi10Header_;
    // bool alphaRestart_;
    // surfaceScalarField alphaPhi10_;
    autoPtr<compressibleInterPhaseTransportModel> turbulence_;
    volScalarField K_;

    fv::options& fvOptions_;
    volScalarField& p_;
    volScalarField& T_;
    const volScalarField& psi1_;
    const volScalarField& psi2_;

    // Alpha fluxes

    IOobject alphaPhi10Header;
    bool alphaRestart;
    surfaceScalarField alphaPhi10;
    surfaceScalarField alphaPhiUn;
    tmp<surfaceScalarField> talphaPhi1Corr0;

    IOMRFZoneList MRF_;



public:

    TypeName("compressibleInterFoam");

    //- Constructors

        compressibleInterFoam
        (
            dynamicFvMesh& mesh_
        );

        //- Destructor
        virtual ~compressibleInterFoam(){}

    //- Member Functions

        //- Solve thermal-hydraulics
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        virtual scalar maxDeltaT();

        //- Correct fluid and structure temperatures without
        //  updating the flow regime map
        void correctEnergy();

        //- Correct pressure, velocities without updating
        //  the flow regime map
        void correctFluidMechanics();

        //- Update and print Courant number
        void correctCourant();

        //- Update continuity errors of the fluid
        void correctContErr();

        //- Print continuity error information
        void printContErr();

        //- Computes the cumulated continuity error
        void calcCumulContErr();

        virtual scalar getResidual()
        {
            return scalar(0);
        };

        virtual void correctBaffleLessFields()
        {}
        virtual void deformMesh()
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicsModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
