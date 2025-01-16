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
    Foam::regionSolvers

Description
    Class to hold the lists of region meshes and solvers

    Also provides loop and iteration functionality which automatically set the
    region Info prefix for each of the solvers returned by the '[]' operator or
    iterator.

Usage
    Given the \c regionSolvers named solvers:
    \verbatim
        // Create the region meshes and solvers
        regionSolvers solvers(runTime);
    \endverbatim

    The list of solvers can be looped over:
    \verbatim
        forAll(solvers, i)
        {
            solvers[i].momentumPredictor();
        }
    \endverbatim
    where the '[]' operator sets the region Info prefix.  After the loop the
    region Info prefix remains set to the last region prefix and so for global
    messages, e.g. the global time-step the Info prefix must be specifically
    reset to spaces by calling the \c setGlobalPrefix() function.

    Alternatively the list of solvers can be iterated over:
    \verbatim
        forAllIter(regionSolvers, solvers, solver)
        {
            solver->momentumPredictor();
        }
    \endverbatim
    where the iterator increment sets the region \c Info prefix and at the end
    automatically resets the \c Info prefix to spaces.

    This class creates the solvers at the highest level. If a solver is a multi-
    physics solver, more solver will be created iteratively. When setting up 
    the simulation, in controlDict the regions are specified as:
    \verbatim
        regionSolvers
        {
            regionName      solverName;
            fluidRegion     onePhase;
            neutroRegion    diffusionNeutronics;
            loop1           multiPhysicsSolver;
        }
    \endverbatim

SourceFiles
    regionSolvers.C

\*---------------------------------------------------------------------------*/

#ifndef regionSolvers_H
#define regionSolvers_H

#include "solver.H"
#include "meshHandler.H"
#include "meshToMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionSolvers Declaration
\*---------------------------------------------------------------------------*/

class regionSolvers
{
    // Private Member Data

        //- Mapper
        meshHandler meshHandler_;

        //- Region solvers
        PtrList<solver> solvers_;

        //- Global space padding prefix
        string prefix0_;

        //- List of space padded region prefixes
        stringList prefixes_;

        bool multiRegionTightCoupled_;
        scalar multiRegionResidual_;
        scalar maxMultiRegionIterations_;

        SquareMatrix<scalar> invRBFmatrix_;


public:

    //Outer loop controls


    // Constructors

        //- Construct from components
        regionSolvers(const Time& runTime);

        //- Disallow default bitwise copy construction
        regionSolvers(const regionSolvers&) = delete;


    //- Destructor
    ~regionSolvers();


    // Member Functions

        bool tightlyCoupled(){return multiRegionTightCoupled_;};
        scalar multiRegionResidual(){return multiRegionResidual_;};
        scalar iterMax() {return maxMultiRegionIterations_;};

        //- Return the number of region solvers
        inline label size() const;

        //- Set the Info prefix to space padding for global messages
        void setGlobalPrefix() const;

        //- Set the Info prefix to the space padded region name
        void setPrefix(const label i) const;

        //- Reset the Info prefix to null
        void resetPrefix() const;

        //void mapFields(const Time& runTime);

        //void initializeMappedFields(const Time& runTime);

        //void createBaffleLessMeshes(const Time& runTime);

        //void createCouplingFields(const Time& runTime);

        meshHandler& mapper()
        {
            return meshHandler_;
        }

        //void interpolateAndMapFields(const Time& runTime);


    // Iterator

        class iterator;
        friend class iterator;

        // Iterator class to loop over the solvers, setting the prefix for each
        class iterator
        {
            // Private Data

                //- Reference to the regionSolvers
                regionSolvers& regionSolvers_;

                //- Current solver index
                label index_;

        public:

            // //Access to meshes

            friend class regionSolvers;

            // Constructors

                //- Construct for the regionSolvers
                inline explicit iterator(regionSolvers&);

            // Member Operators

                inline bool operator==(const iterator&) const;
                inline bool operator!=(const iterator&) const;

                inline solver& operator*();
                inline solver& operator()();
                inline solver* operator->() const;

                inline iterator operator++();
                inline iterator operator++(int);
        };

        //- Return an iterator to begin traversing the solvers
        inline iterator begin();

        //- Return an iterator to end traversing the solvers
        inline iterator end();


    // Member Operators

        //- Cast to the list of solvers
        inline operator PtrList<solver>&();

        //- Set the region i prefix and return the corresponding solver
        solver& operator[](const label i);

        //- Disallow default bitwise assignment
        void operator=(const regionSolvers&) = delete;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regionSolversI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
