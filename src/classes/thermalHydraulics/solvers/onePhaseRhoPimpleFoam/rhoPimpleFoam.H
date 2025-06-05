/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
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
    Foam::thermalHydraulicsModels::rhoPimpleFoam

Description
    This class encapsulates a single fluid + 1 stationary structure
    Eulerian solver with a user-defineable regime-map based phyiscs modelling.

SourceFiles
    rhoPimpleFoam.C

\mainauthor
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

\*---------------------------------------------------------------------------*/

#ifndef rhoPimpleFoam_H
#define rhoPimpleFoam_H

#include "dynamicFvMesh.H"
#include "turbulentFluidThermoModel.H"
#include "CorrectPhi.H"
#include "pressureControl.H"
#include "solver.H"
#include "IOMRFZoneList.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "volPointInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
//class regime;
//class regimeMapModel;

namespace solvers
{

/*---------------------------------------------------------------------------*\
                      Class rhoPimpleFoam Declaration
\*---------------------------------------------------------------------------*/

class rhoPimpleFoam
:     
    public solver
{
private:


    dynamicFvMesh& mesh_;
    pimpleControl pimple_;
    autoPtr<fluidThermo> pThermo_;
    fluidThermo& thermo_;

    volScalarField p_;
    volScalarField rho_;
    volVectorField U_;
    surfaceScalarField phi_;
    
    pressureControl pressureControl_;
    autoPtr<compressible::turbulenceModel> turbulence_;

    volScalarField dpdt_;
    volScalarField K_;

    IOMRFZoneList MRF_;
    const dimensionedScalar rhoMax_;
    const dimensionedScalar rhoMin_;

    fv::options& fvOptions_;

    bool correctPhi_;
    bool moveMeshOuterCorrectors_;
    bool checkMeshCourantNo_;

    autoPtr<surfaceVectorField> rhoUf_;

    scalar residual_;
    scalar cumulativeContErr_;

    //- Original mesh
    pointField originalPoints_;
    



public:

    TypeName("rhoPimpleFoam");

    //- Constructors

        rhoPimpleFoam
        (
            dynamicFvMesh& mesh
        );

        //- Destructor
        virtual ~rhoPimpleFoam(){}

    //- Member Functions

        //- Solve thermal-hydraulics
        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        virtual scalar maxDeltaT();

        //- Update continuity errors of the fluid
        void compressibleContinuityErrs();

        virtual scalar getResidual()
        {
            return residual_;
        };

        virtual void correctBaffleLessFields();

        template<class Type>
        void correctBaffleLessField(word fieldName, fvMesh& baffleLessMesh);

        void deformMesh();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace thermalHydraulicsModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
