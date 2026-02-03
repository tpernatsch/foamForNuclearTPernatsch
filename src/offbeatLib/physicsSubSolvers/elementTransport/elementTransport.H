/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::elementTransport

Description
    Mother class for handling element transport solvers.

    \htmlinclude elementTransport_elementTransport.html

    \warning
    For backward-compatibility, `fromLatestTime` can be still used as the typename
    for the `none` elementTransport model. However, its use is deprecated.

SourceFiles
    elementTransport.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef elementTransport_H
#define elementTransport_H

#include "dynamicFvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"
#include "primitiveFields.H"
#include "materials.H"
#include "transportSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class elementTransport Declaration
\*---------------------------------------------------------------------------*/

//class materialModel;

class elementTransport
{
    // Private data

    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        elementTransport(const elementTransport&);

        //- Disallow default bitwise assignment
        void operator=(const elementTransport&);

protected:
    
    // Protected data

        //- Reference to the mesh
        fvMesh& mesh_;

        //- Const reference to the materials class
        const materials& mat_;

        //- Element transport dictionary
        const dictionary elementTransportDict_;
        
    // Protected member functions

public:

    //- Runtime type information
    TypeName("none");
    

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            elementTransport,
            dictionary,
            (
                fvMesh& mesh, 
                const materials& mat,
                const dictionary& elementTransportDict
            ),
            (mesh, mat, elementTransportDict)
        );


    // Constructors

        //- Construct from material dictionary
        elementTransport
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict
        );

        //- Return a pointer to a new elementTransport
        static autoPtr<elementTransport> New
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& solverDict
        );


    //- Destructor
    virtual ~elementTransport();


    // Access Functions

    // Member Functions

        //- Correct
        virtual void correct(){};

        //- Converged
        virtual bool converged()
        {
            return true;
        };

        //- 
        virtual scalar nextDeltaT()
        {
            return GREAT;
        };

        virtual void updateMesh(){};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
