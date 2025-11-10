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
    Foam::loop

Description
    This class creates a multiphysics tightly coupled loop, made of an arbitrary
    number of physics. The class member-functions are inherited from solver.H
    and are the same as those of the other solvers. Each member function consists
    of a loop where the same member function is called iteratively for all the
    solvers included in the loop. Note that another loop
    can be created within a loop, hence creating a hierarchical
    structure which might be helpful for multi-scale problems.

Usage
    In regionsDict, the regionSolvers subDict is used to create a
    loop

    \verbatim
        regionSolvers
        {
            multiPhysicsSolver1 //multiPhysicsSolverName
            {
                solvers // list of solvers
                {
                    regionName              solverName;
                    fluidRegion             onePhase;
                    neutroRegion            diffusionNeutronics;
                    multiPhysicsSolver2     loop;
                }
                minResidual 0.00005;   // residual after which the iteration is ended
                maxIterations 3;       // maximum number of iterations before the loop is ended
                includeFluidMechanicsInLoop  false; //defaulted to false
            }
            multiPhysicsSolver2
            {
                ...
            }
        }

    \endverbatim

    If includeFluidMechanicsInLoop entry is set to false (default value), for the thermal-hydraulics
    solvers only the energy is solved iteratively with the other physics, while the fluid-mechanics
    is solved only once at the beginning of the iteration. This is done so that useless
    fluid-mechanics iterations are solved.

SourceFiles
    loop.C

\mainauthor Author of this file (and associated .C or included .H files):
    Giovanni Nervi <giovanni.nervi@epfl.ch>, EPFL (Switzerland)

\*---------------------------------------------------------------------------*/

#ifndef loop_H
#define loop_H

#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"
#include "solver.H"
#include "meshHandler.H"
#include "dynamicFvMesh.H"
#include "volPointInterpolation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class loop Declaration
\*---------------------------------------------------------------------------*/

class loop
:
    public solver,
    public IOdictionary
{

protected:

    // Protected data

        //- Reference to mesh
        PtrList<solver> solvers_;

        //-
        meshHandler* meshHandler_;

        //-
        wordList solverNames_;

        //-
        wordList singlePhysicsSolverNames_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        loop(const loop&);

        //- Disallow default bitwise assignment
        void operator=(const loop&);


public:

    //- Runtime type information
    TypeName("loossssp")

    // Constructors

        //- Construct from mesh
        loop(dynamicFvMesh& mesh);

    // Selectors

        static autoPtr<loop> New
        (
            dynamicFvMesh& mesh
        );


    //- Destructor
    	virtual ~loop(){}



    // Member Functions

        //- Needed to remove abstractedness from IOdictionary
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }

        //- Access

        virtual void correctBaffleLessFields();

        virtual void deformMesh();

        virtual void correctPhysics() = 0;
        virtual void correctTightlyCoupledPhysics(){};

        virtual scalar maxDeltaT() ;

        virtual scalar getResidual();

        virtual void getMapper(meshHandler& mapper)
        {
            meshHandler_ =&mapper;
        };

        virtual void createSolvers(word name) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
