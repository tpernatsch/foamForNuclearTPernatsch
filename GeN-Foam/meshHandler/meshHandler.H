/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2022-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::meshHandler

Description
    Class to hold the lists of region meshes and solvers

    Also provides loop and iteration functionality which automatically set the
    region Info prefix for each of the solvers returned by the '[]' operator or
    iterator.

Usage
    Given the \c meshHandler named solvers:
    \verbatim
        // Create the region meshes and solvers
        meshHandler solvers(runTime);
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
        forAllIter(meshHandler, solvers, solver)
        {
            solver->momentumPredictor();
        }
    \endverbatim
    where the iterator increment sets the region \c Info prefix and at the end
    automatically resets the \c Info prefix to spaces.

SourceFiles
    meshHandler.C

\*---------------------------------------------------------------------------*/

#ifndef mapping_H
#define mapping_H

#include "solver.H"
#include "meshToMesh.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class meshHandler Declaration
\*---------------------------------------------------------------------------*/

class meshHandler
{
    // Private Member Data

        //- Region meshes
        
        PtrList<fvMesh> meshes_;

        PtrList<fvMesh> mappingMeshes_;

        dictionary mappingDict_;
        PtrList<PtrList<meshToMesh>> mappingList_;
        PtrList<PtrList<volScalarField>> scalarCouplingFields_;
        PtrList<PtrList<volVectorField>> vectorCouplingFields_;


public:

    //Outer loop controls


    // Constructors

        //- Construct from components
        meshHandler(const Time& runTime);

        //- Disallow default bitwise copy construction
        meshHandler(const meshHandler&) = delete;


    //- Destructor
    ~meshHandler();


    // Member Functions

        void mapAllFields(const Time& runTime);

        void mapTheseFields(const Time& runTime, wordList meshToMap);

        void initializeMappedFields(const Time& runTime);

        void createBaffleLessMeshes(const Time& runTime);

        void createCouplingFields(const Time& runTime);

        fvMesh& returnMesh(word meshName);
        fvMesh& returnMappingMesh(word meshName);
        bool contains(wordList list, word thisWord);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
