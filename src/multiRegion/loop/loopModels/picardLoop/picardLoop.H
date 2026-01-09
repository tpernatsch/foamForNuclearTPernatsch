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
    Foam::picardLoop

Description
    This class creates a multiphysics tightly coupled loop, made of an arbitrary
    number of physics. The class member-functions are inherited from solver.H
    and are the same as those of the other solvers. Each member function consists
    of a loop where the same member function is called iteratively for all the
    solvers included in the picardLoop. Note that another picardLoop
    can be created within a picardLoop, hence creating a hierarchical
    structure which might be helpful for multi-scale problems.

\par Options

    \table
        Property        | Description               | Req'd | Default value
        subSolvers      | List of subsolvers        | yes   |
        maxIterations   | Maximum number of iterations | yes |
        minResidual     | Maximum residual          | yes |
        couplingStartTime | Time after which the coupling starts | no | 0
    \endtable

Usage
    In regionsDict, the following subDict is used to create a picardLoop

    \verbatim
        regionSolvers
        {
            Level_0
            {
                Level_1     picardLoop;
            }
            Level_1
            {
                subSolvers // list of solvers
                {
                    regionName              solverName;
                    fluidRegion             onePhase;
                    neutroRegion            diffusionNeutronics;
                    Level_2                 picardLoop;
                }
                minResidual 0.00005;   // residual after which the iteration is ended
                maxIterations 3;       // maximum number of iterations before the loop is ended
            }

            Level_2
            {
                ...
            }
        }
    \endverbatim

    The picardLoop automatically solves all the equations of the solvers for each time step.
    In some cases, it is convenient to only solve the energy equation of the
    thermal-hydraulics solver, while the fluid-mechanics is solved only once at the
    beginning of the iteration. This can be done by setting the picardLoopNoFluid solover instead.

SourceFiles
    picardLoop.C

\mainauthor Author of this file (and associated .C or included .H files):
    Giovanni Nervi <giovanni.nervi@epfl.ch>, EPFL (Switzerland)

\*---------------------------------------------------------------------------*/

#ifndef picardLoop_H
#define picardLoop_H

#include "volFields.H"
#include "fvCFD.H"
#include "meshHandler.H"
#include "dynamicFvMesh.H"
#include "volPointInterpolation.H"
#include "solver.H"
#include "loop.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{
namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class picardLoop Declaration
\*---------------------------------------------------------------------------*/

class picardLoop
:
    public loop
{

protected:

    // Protected data

        //-
        dictionary multiPhysicsDict_;

        //-
        dictionary loopProperties_;

        //-
        scalar maxResidual_;

        //-
        label maxIterations_;

        //-
        scalar couplingStartTime_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        picardLoop(const picardLoop&);

        //- Disallow default bitwise assignment
        void operator=(const picardLoop&);


public:

    //- Runtime type information
    TypeName("picardLoop")


    // Constructors

        //- Construct from mesh
        picardLoop(dynamicFvMesh& mesh);

    // Selectors


    //- Destructor
    	virtual ~picardLoop(){}


    // Member Functions

        // Access

        //-
        virtual void createSolvers(word name);

        //-
        virtual void correctPhysics();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
