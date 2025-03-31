/*---------------------------------------------------------------------------*\
|                                                              |
|      / /      / | / /          / /               |
|     / /    /  \  /  |/ /    / /     /  \ /  `/  /  ` \    |
|    / // /  /  / / /|  /  // / /    / // // // /  / / / / / /    |
|    \/   \/ // |/          //       \/ \,/  // // //     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2412                                                  |
|    Copyright 2011-2016 OpenFOAM Foundation, 2017-2024 OpenCFD Ltd.          |
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
    Foam::thermalHydraulicsModels::compressibleInterDyMFoam

Description
    Solver for two compressible, non-isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach. This
    standard solver has been transposed from OpenFOAM.

SourceFiles
    compressibleInterDyMFoam.C

\*---------------------------------------------------------------------------*/

#ifndef compressibleInterDyMFoamH
#define compressibleInterDyMFoamH

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
#include "dynamicFvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
//class regime;
//class regimeMapModel;

namespace solvers
{

/*---------------------------------------------------------------------------*\
                      Class compressibleInterDyMFoam Declaration
\*---------------------------------------------------------------------------*/

class compressibleInterDyMFoam
:
    public solver
{
private:

    //- Mesh
    dynamicFvMesh& dynMesh;

    //-
    pimpleControl pimple;

    //- Phases

    scalar cumulativeContErr;


    volScalarField p_rgh;
    volVectorField U;
    surfaceScalarField phi;
    twoPhaseMixtureThermo mixture;
    volScalarField& alpha1;
    volScalarField& alpha2;
    const volScalarField& rho1;
    const volScalarField& rho2;
    volScalarField rho;
    dimensionedScalar pMin;

    //- Gravity related fields
    uniformDimensionedVectorField g;
    uniformDimensionedScalarField hRef;
    dimensionedScalar ghRef;
    volScalarField gh;
    surfaceScalarField ghf;

    surfaceScalarField rhoPhi;
    volScalarField dgdt;

    // //Create alphaFluxes
    // IOobject alphaPhi10Header;
    // bool alphaRestart;
    // surfaceScalarField alphaPhi10;
    autoPtr<compressibleInterPhaseTransportModel> turbulence;
    surfaceVectorField Uf;
    volScalarField K;

    fv::options& fvOptions;
    volScalarField& p;
    volScalarField& T;
    const volScalarField& psi1;
    const volScalarField& psi2;

    // Alpha fluxes

    IOobject alphaPhi10Header;
    bool alphaRestart;
    surfaceScalarField alphaPhi10;
    surfaceScalarField alphaPhiUn;
    tmp<surfaceScalarField> talphaPhi1Corr0;

    IOMRFZoneList MRF;



public:

    TypeName("compressibleInterDyMFoam");

    //- Constructors

        compressibleInterDyMFoam
        (
            dynamicFvMesh& mesh
        );

        //- Destructor
        virtual ~compressibleInterDyMFoam(){}

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
