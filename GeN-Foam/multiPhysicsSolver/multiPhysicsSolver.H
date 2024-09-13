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
    Base class for thermo-mechanical models.

    Reference publications:

    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM� based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226

    Author of this file (and associated .C or included .H files): 
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>, EPFL (Switzerland)

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
        meshHandler* meshHandler_;
        dictionary multiPhysicsDict_;
        wordList solverNames_;
        scalar minResidual_;
        label maxIterations_;


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
