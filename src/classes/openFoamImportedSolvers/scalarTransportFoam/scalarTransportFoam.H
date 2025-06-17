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
    Foam::thermalHydraulicsModels::scalarTransportFoam

Description
    OpenFOAM-based passive scalar transport equation solver
    ([link to OpenFOAM documentation](https://www.openfoam.com/documentation/guides/latest/doc/guide-applications-solvers-basic-scalarTransportFoam.html)).

    The equation is given by:

    $$
        \frac{\partial T}{\partial t} + \nabla \cdot \left( \mathbf{u} T \right) - \nabla \cdot \left(D_T \nabla T \right)
        = S_T
    $$

    \vartable
    T       | Passive scalar
    D_T     | Diffusion coefficient
    S_T     | Source
    \endvartable

    \heading Required fields
    \plaintable
        T       | Passive scalar
        U       | Velocity [m/s]
    \endplaintable


SourceFiles
    scalarTransportFoam.C

\*---------------------------------------------------------------------------*/

#ifndef scalarTransportFoam_H
#define scalarTransportFoam_H

#include "fvCFD.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "solver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                      Class scalarTransportFoam Declaration
\*---------------------------------------------------------------------------*/

class scalarTransportFoam
:
    public solver
{
private:

    //- Mesh
    fvMesh& mesh_;

    //-
    pimpleControl pimple_;

    //- Scalar field to be transported
    volScalarField T_;

    //- Velocity field
    volVectorField U_;

    //-
    surfaceScalarField phi_;

    //-
    IOdictionary transportProperties_;

    //-
    dimensionedScalar DT_;

    //-
    fv::options& fvOptions_;


public:

    TypeName("scalarTransportFoam");

    //- Constructors

        scalarTransportFoam
        (
            fvMesh& mesh_
        );

        //- Destructor
        virtual ~scalarTransportFoam(){}

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
