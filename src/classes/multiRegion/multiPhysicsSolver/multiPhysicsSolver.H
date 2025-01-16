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
    Foam::multiPhysicsSolver

Description
    This class creates a multiphysics tightly coupled loop, made of an arbitrary
    number of physics. The class member-functions are inherited from solver.H 
    and are the same as those of the other solvers. Each member function consists
    of a loop where the same member function is called iteratively for all the 
    solvers included in the multiPhysicsSolver. Note that another multiphysicsSolver
    can be created within a multiphysicsSolver, hence creating a hierarchical
    structure which might be helpful for multi-scale problems.

Usage
    In multiRegionCouplingDict, the multiPhysicsSolvers subDict is used to create a
    multiPhysicsSolver

    \verbatim
        multiPhysicsSolvers
        {
            multiPhysicsSolver1 //multiPhysicsSolverName
            {
                solvers // list of solvers
                {
                    regionName              solverName;
                    fluidRegion             onePhase;
                    neutroRegion            diffusionNeutronics;
                    multiPhysicsSolver2     multiPhysicsSolver;
                }
                minResidual 0.00005;   // residual after which the iteration is ended
                maxIterations 3;       // maximum number of iterations before the loop is ended
                includeFluidMechanicsInLoop  false; //defaulted to false
            }
        }

        multiPhysicsSolver2
        {
            ...   
        }
    \endverbatim
    }

    If includeFluidMechanicsInLoop entry is set to false (default value), for the thermal-hydraulics
    solvers only the energy is solved iteratively with the other physics, while the fluid-mechanics
    is solved only once at the beginning of the iteration. This is done so that useless 
    fluid-mechanics iterations are solved.

SourceFiles
    multiPhysicsSolver.C






\*---------------------------------------------------------------------------*/

#ifndef thermoMechanics_H
#define thermoMechanics_H

#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"
#include "solver.H"
#include "meshHandler.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace solvers
{

/*---------------------------------------------------------------------------*\
                          Class multiPhysicsSolver Declaration
\*---------------------------------------------------------------------------*/

class multiPhysicsSolver
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
        dictionary multiPhysicsDict_;

        //-
        wordList solverNames_;

        //-
        wordList singlePhysicsSolverNames_;

        //-
        scalar minResidual_;

        //-
        label maxIterations_;

        //- 
        List<bool> correctOnlyEnergy_; 


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        multiPhysicsSolver(const multiPhysicsSolver&);

        //- Disallow default bitwise assignment
        void operator=(const multiPhysicsSolver&);


public:

    //- Runtime type information
    TypeName("multiPhysicsSolver")

// #ifndef SWIG
//         declareRunTimeSelectionTable
//         (
//             autoPtr,
//             multiPhysicsSolver,
//             dictionary,
//             (
//                 fvMesh& mesh,
//                 bool multiRegion
//             ),
//             (mesh)
//         );
// #endif

    // Constructors

        //- Construct from mesh
        multiPhysicsSolver(fvMesh& mesh);


    //- Destructor
    	virtual ~multiPhysicsSolver(){}


    // Member Functions

        //- Needed to remove abstractedness from IOdictionary
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }

        //- Access

        virtual void correctBaffleLessFields();

        virtual void deformMesh();

        virtual void correctPhysics();
        virtual void correctTightlyCoupledPhysics();

        virtual scalar maxDeltaT();

        virtual scalar getResidual();

        virtual void getMapper(meshHandler& mapper)
        {
            meshHandler_ =&mapper;
        };

        void createSolvers(word name);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace solvers
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
