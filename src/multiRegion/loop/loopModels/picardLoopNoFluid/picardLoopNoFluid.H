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
    Foam::picardLoopNoFluid

Description
    Identical to the picardLoop, but it does not solve the fluid-mechancs equations at
    each of the iterations. Instead, it solves the fluid-mechanics equations only once
    at the beginning of the iteration, and then solves the energy equation of the
    thermal-hydraulics solver iteratively with the other physics.

\par Options

    \table
        Property        | Description               | Req'd | Default value
        subSolvers      | List of subsolvers        | yes   |
        maxIterations   | Maximum number of iterations | yes |
        maxResidual     | Maximum residual          | yes |
        couplingStartTime | Time after which the coupling starts | no | 0
    \endtable

Usage
    In regionsDict, the following subDict is used to create a picardLoopNoFluid

    \verbatim
        regionSolvers
        {
            Level_0
            {
                Level_1     picardLoopNoFluid;
            }
            Level_1
            {
                subSolvers // list of solvers
                {
                    regionName              solverName;
                    fluidRegion             onePhase;
                    neutroRegion            diffusionNeutronics;
                    Level_2                 picardLoopNoFluid;
                }
                maxResidual 0.00005;   // residual after which the iteration is ended
                maxIterations 3;       // maximum number of iterations before the loop is ended
            }

            Level_2
            {
                ...
            }
        }
    \endverbatim

SourceFiles
    picardLoopNoFluid.C

\mainauthor Author of this file (and associated .C or included .H files):
    Giovanni Nervi <giovanni.nervi@epfl.ch>, EPFL (Switzerland)

\*---------------------------------------------------------------------------*/

#ifndef picardLoopNoFluid_H
#define picardLoopNoFluid_H

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
                          Class picardLoopNoFluid Declaration
\*---------------------------------------------------------------------------*/

class picardLoopNoFluid
:
    public loop
{

protected:

    // Protected data

        //-

        dictionary multiPhysicsDict_;
        dictionary loopProperties_;

        scalar maxResidual_;
        label maxIterations_;
        scalar couplingStartTime_;




private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        picardLoopNoFluid(const picardLoopNoFluid&);

        //- Disallow default bitwise assignment
        void operator=(const picardLoopNoFluid&);


public:

    //- Runtime type information
    TypeName("picardLoopNoFluid")


    // Constructors

        //- Construct from mesh
        picardLoopNoFluid(dynamicFvMesh& mesh);

    // Selectors


    //- Destructor
    	virtual ~picardLoopNoFluid(){}


    // Member Functions

        //- Access

        virtual void createSolvers(word name);
        virtual void correctPhysics();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
