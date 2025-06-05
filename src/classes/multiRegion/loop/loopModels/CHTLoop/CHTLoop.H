/*---------------------------------------------------------------------------*\
|       ______          _   __           ______                               |
|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |
|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |
|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |
|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |
|    Copyright (C) 2015 - 2022 EPFL                                           |
|                                                                             |
|    Built on OpenFOAM v2212                                                  |
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
    Foam::CHTLoop

Description
    This class creates a multiphysics tightly coupled loop, made of an arbitrary
    number of physics. The class member-functions are inherited from solver.H
    and are the same as those of the other solvers. Each member function consists
    of a loop where the same member function is called iteratively for all the
    solvers included in the CHTLoop. Note that another CHTLoop
    can be created within a CHTLoop, hence creating a hierarchical
    structure which might be helpful for multi-scale problems.

Usage
    In multiRegionCouplingDict, the CHTLoops subDict is used to create a
    CHTLoop

    \verbatim
        CHTLoops
        {
            CHTLoop1 //CHTLoopName
            {
                solvers // list of solvers
                {
                    regionName              solverName;
                    fluidRegion             onePhase;
                    neutroRegion            diffusionNeutronics;
                    CHTLoop2     CHTLoop;
                }
                minResidual 0.00005;   // residual after which the iteration is ended
                maxIterations 3;       // maximum number of iterations before the loop is ended
                includeFluidMechanicsInLoop  false; //defaulted to false
            }
        }

        CHTLoop2
        {
            ...
        }
    \endverbatim

    If includeFluidMechanicsInLoop entry is set to false (default value), for the thermal-hydraulics
    solvers only the energy is solved iteratively with the other physics, while the fluid-mechanics
    is solved only once at the beginning of the iteration. This is done so that useless
    fluid-mechanics iterations are solved.

SourceFiles
    CHTLoop.C

\*---------------------------------------------------------------------------*/

#ifndef CHTLoop_H
#define CHTLoop_H

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
                          Class CHTLoop Declaration
\*---------------------------------------------------------------------------*/

class CHTLoop
:
    public loop
{
protected:

    // Protected data

        //-

        dictionary multiPhysicsDict_;
        dictionary CHTProperties_; 

        word fluidRegionName_;
        label fluidPatchID_;
        word solidRegionName_;
        label solidPatchID_;
        label fluidSolverID_;
        label solidSolverID_;
        scalar couplingStartTime_;
        scalar minResidual_;
        label maxIterations_;
        bool useHTC_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        CHTLoop(const CHTLoop&);

        //- Disallow default bitwise assignment
        void operator=(const CHTLoop&);


public:

    //- Runtime type information
    TypeName("CHTLoop")


    // Constructors

        //- Construct from mesh
        CHTLoop(dynamicFvMesh& mesh);

    // Selectors


    //- Destructor
    	virtual ~CHTLoop(){}


    // Member Functions

        //- Access
        virtual void createSolvers(word name);
        virtual void correctPhysics();

        void fromSolidToFluid();
        void fromFluidToSolid();
        scalar calcFSIResidual();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
