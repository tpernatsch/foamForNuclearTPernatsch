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
    Foam::FSILoop

Description
    This class implements a dedicated Fluid–Structure Interaction (FSI) solver
    that couples a single fluid region to a single solid region using the
    Arbitrary Lagrangian–Eulerian (ALE) approach. The solver is designed as a
    specialized loop that advances both physics in a tightly coupled manner
    within each time step, ensuring consistent interface conditions for velocity,
    displacement, and traction. It assumkes that the solid and fluid are coupled
    through a single patch on both sides, allowing for a direct exchange of
    information between the fluid and solid regions. The fluid solver is expected 
    to support dynamicFvMesh in order to deform the mesh according to the solid 
    displacement. The loop assumes that the solid region is not deformed, using 
    a non-updated approach for the displacement solution.

Usage
    In regionsDict, a subDict in regionsDict/regionSolvers is used to create a FSILoop

    \verbatim
        Level_1
        {
            subSolvers
            {
                fluid rhoPimpleFoam;
                solid extendedThermoMechanics;
                ...
                arbitraryName  solverName;
            }
            
            maxIterations 100;
            minResidual 1e-5;
            
            FSILoopCoeffs
            {
                fluidRegionName fluid; // Name of the fluid region
                solidRegionName solid; // Name of the solid region
                solidSidePatchName plate; // Name of the coupled patch on the solid side
                fluidSidePatchName plate; // Name of the coupled patch on the fluid side
                underRelaxationFactor 0.05; // Under-relaxation factor for the coupling
                couplingStartTime 2; // Time after which the coupling starts
            }
        }
    \endverbatim

SourceFiles
    FSILoop.C

\*---------------------------------------------------------------------------*/

#ifndef FSILoop_H
#define FSILoop_H

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
                          Class FSILoop Declaration
\*---------------------------------------------------------------------------*/

class FSILoop
:
    public loop
{
protected:

    // Protected data

        //-

        dictionary multiPhysicsDict_;
        dictionary FSIProperties_; 

        word fluidRegionName_;
        label fluidPatchID_;
        word solidRegionName_;
        label solidPatchID_;
        label fluidSolverID_;
        label solidSolverID_;
        scalar underRelaxation_;
        scalar couplingStartTime_;
        scalar minResidual_;
        label maxIterations_;
        bool thermalCoupling_;
        
        autoPtr<volPointInterpolation> solidInterpolator_;
        autoPtr<vectorField> oldDisplacementAtFaces_;
        autoPtr<volVectorField> fluidDisplacementAtSolid_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        FSILoop(const FSILoop&);

        //- Disallow default bitwise assignment
        void operator=(const FSILoop&);


public:

    //- Runtime type information
    TypeName("FSILoop")


    // Constructors

        //- Construct from mesh
        FSILoop(dynamicFvMesh& mesh);

    // Selectors


    //- Destructor
    	virtual ~FSILoop(){}


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
