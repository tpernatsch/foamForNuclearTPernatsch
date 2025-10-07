/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2506                                                  |
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
    Foam::solver

Description
    Abstract base class for run-time selectable region solvers.

SourceFiles
    solver.C
    solverNew.C

\*---------------------------------------------------------------------------*/

#ifndef solver_H
#define solver_H

#include "dynamicFvMesh.H"
#include "argList.H"
// #include "pimpleNoLoopControl.H"
// #include "fvModels.H"
// #include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solver Declaration
\*---------------------------------------------------------------------------*/

class solver
:
    public regIOobject
{
protected:

    // Protected Data

        //- Region mesh
        dynamicFvMesh& mesh_;

        //- Switch for steady-state
        //  Set automatically from the selected ddt scheme
        bool steady;

        //- Switch for local time step transient operation
        //  Set automatically from the selected ddt scheme
        bool LTS;


private:

    // Private Member Data

        //- Pointer to the demand driven fvModels MeshObject
        // mutable Foam::fvModels* fvModelsPtr;

        //- Pointer to the demand driven fvConstraints MeshObject
        // mutable Foam::fvConstraints* fvConstraintsPtr;


    // Private Member Functions

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream&) const;


public:

    // Public Data

        //- Region mesh
        const dynamicFvMesh& mesh;

        //- Time
        const Time& runTime;

        //- PIMPLE inner-loop controls
        // pimpleNoLoopControl pimple;

        //- deltaT increase factor
        static scalar deltaTFactor;


    //- Runtime type information
    TypeName("solver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solver,
            dynamicFvMesh,
            (dynamicFvMesh& mesh),
            (mesh)
        );


    // Constructors

        //- Construct from region mesh
        solver(dynamicFvMesh& mesh);

        //- Disallow default bitwise copy construction
        solver(const solver&) = delete;


    // Selectors

        //- Select, construct and return the solver
        static autoPtr<solver> New
        (
            const word& solverName,
            dynamicFvMesh& mesh
        );

        //- Load the specified solver library
        static void load(const word& solverName);


    //- Destructor
    virtual ~solver();


    // Member Functions

        bool transient() const
        {
            return !steady && !LTS;
        }

        //- Return the fvModels that are created on demand
        // Foam::fvModels& fvModels() const;

        // //- Return the fvConstraints that are created on demand
        // Foam::fvConstraints& fvConstraints() const;

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() = 0;

        //- Called at the start of the time-step, before the PIMPLE loop
        // virtual void preSolve() = 0;

        // //- Called at the start of the PIMPLE loop to move the mesh
        // virtual void moveMesh() = 0;

        // //- Called at the start of the PIMPLE loop
        // virtual void prePredictor() = 0;

        // //- Construct and optionally solve the momentum equation
        // virtual void momentumPredictor() = 0;

        // //- Construct and solve the thermophysical property equations,
        // //  convert to temperature
        // //  and update thermophysical and transport properties
        // virtual void thermophysicalPredictor() = 0;

        // //- Construct and solve the pressure equation in the PISO loop
        // virtual void pressureCorrector() = 0;

        // //- Correct the momentum and thermophysical transport modelling
        // virtual void postCorrector() = 0;

        // //- Called after the PIMPLE loop at the end of the time-step
        // virtual void postSolve() = 0;

        virtual void correctPhysics() = 0;
        virtual void correctTightlyCoupledPhysics() = 0;

        virtual scalar getResidual()=0;

        //-BaffleLess mesh handling
        virtual void correctBaffleLessFields()=0;

        virtual void deformMesh()=0;

        word regionName()
        {
            return mesh_.name();
        }


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const solver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
