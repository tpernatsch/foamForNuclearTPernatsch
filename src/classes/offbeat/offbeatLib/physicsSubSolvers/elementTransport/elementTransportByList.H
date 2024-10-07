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
    Foam::elementTransportByList

Description
    Class that gathers all the user-selected transport solvers in a ptrList, 
    to handle the solution of the respective equations.

    \htmlinclude elementTransport_elementTransportByList.html

SourceFiles
    elementTransportByList.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef elementTransportByList_H
#define elementTransportByList_H

#include "elementTransport.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class elementTransportByList Declaration
\*---------------------------------------------------------------------------*/

class elementTransportByList
:
    public elementTransport
{
    // Private data
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        elementTransportByList(const elementTransportByList&);

        //- Disallow default bitwise assignment
        void operator=(const elementTransportByList&);

protected:

        //- PtrList of transportSolver
        PtrList<transportSolver> solverList_;

public:

    //- Runtime type information
    TypeName("byList");
    

    // Constructors

        //- Construct from mesh, materials and dictionary
        elementTransportByList
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict
        ); 
    

    //- Destructor
    ~elementTransportByList();

    // Access Functions

        //- Return reference to constitutive laws list
        inline const PtrList<transportSolver>& solverList() const
        {
            return solverList_;
        }
        
    // Member Functions

        //- Correct
        virtual void correct();

        //- Converged
        virtual bool converged();

        //- Return next delta t for time-stepping criteria
        virtual scalar nextDeltaT();
        
        virtual void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
