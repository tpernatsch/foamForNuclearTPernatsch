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
private:

    // Private Member Data


    // Private Member Functions

        //- Dummy write for regIOobject
        virtual bool writeData(Ostream&) const;
        

protected:

    // Protected Data

        //- Region mesh
        dynamicFvMesh& mesh_;

        //- Time
        const Time& runTime_;

        //- Switch for steady-state
        //  Set automatically from the selected ddt scheme
        bool steady_;

        //- Switch for local time step transient operation
        //  Set automatically from the selected ddt scheme
        //  TODO: what is it? why is it here?
        bool LTS_;


public:

    // Public Data

        //- deltaT increase factor
        static scalar deltaTFactor_;


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

        //- What is it? Why in mother class?
        bool transient() const
        {
            return !steady_ && !LTS_;
        }

        //- Return the current maximum time-step for stable solution
        virtual scalar maxDeltaT() = 0;
        
        //- Return end time
        virtual scalar endTime()
        {
            return GREAT;
        }

        //- Advance physics
        virtual void correctPhysics() = 0;

        //- TODO: Add right comment here
        virtual void correctTightlyCoupledPhysics()
        {
            this->correctPhysics();
        }

        //- BaffleLess mesh handling
        virtual void correctBaffleLessFields(){};

        //- Deform region mesh
        virtual void deformMesh(){};


    // Access Functions

        //- TODO: Why use "get" instead of simply "residual()"
        //  TODO: Add comment
        virtual scalar getResidual()=0;

        //- TODO: Why use "get" instead of simply "residual()"
        //  TODO: Add comment
        virtual scalar getThermalResidual()
        {
            return (getResidual());
        }

        //- Name of region
        word regionName()
        {
            return mesh_.name();
        }

        //- Ref to runTime
        inline virtual const Time& runTime()
        {
            return runTime_;
        }

        //- Ref to mesh
        inline virtual const fvMesh& mesh()
        {
            return mesh_;
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
