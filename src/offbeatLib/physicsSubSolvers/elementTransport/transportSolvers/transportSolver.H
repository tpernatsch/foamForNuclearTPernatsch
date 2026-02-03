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
    Foam::transportSolver

Description
    Mother class for transport solvers.

    \htmlinclude transportSolver_transportSolver.html

SourceFiles
    transportSolver.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, C. Fiorina - EPFL

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef transportSolver_H
#define transportSolver_H

#include "fvMesh.H"
#include "materials.H"
#include "IOdictionary.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class transportSolver Declaration
\*---------------------------------------------------------------------------*/

class transportSolver
{
    // Private data     
        
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        transportSolver(const transportSolver&);

        //- Disallow default bitwise assignment
        void operator=(const transportSolver&);
        
protected:

    // Protected data
        
        //- Reference to mesh
        fvMesh& mesh_;

        //- Const reference to the materials class
        const materials& mat_;

        //- Const reference to the elementTransportDictionary
        const dictionary& dict_;

        //- Residual of latest inner iteration
        scalar residual_;

        //- Residual of initial inner iteration 
        scalar initialResidual_;  

        //- A different residual
        scalar relResidual_;

    // Protected Member Functions    

        //- Read and return number of inner correctors
        int nCorrectors(const word);
        
public: 

    //- Runtime type information
    TypeName("transportSolver");  
    
    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        transportSolver,
        dictionary,
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        ),
        (mesh, mat, elementTransportDict, solverName)
    );
        
    // Constructors

        //- Construct from mesh, materials and dictionary
        transportSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        ); 

        //- Return a pointer to a new transport solver
        static autoPtr<transportSolver> New
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        );

    //- Destructor
    virtual ~transportSolver();

    // Access functions

        //- Return name of the solving variable
        virtual const word& fieldName() = 0;

    // Member Functions 

        //- Correct 
        virtual void correct(){};

        //- Return true if converged
        virtual bool converged( const word& fieldName );

        //- Return true if converged
        virtual bool converged()
        {
            return converged(fieldName());
        };

        //- Next delta t for time stepping criterion
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
