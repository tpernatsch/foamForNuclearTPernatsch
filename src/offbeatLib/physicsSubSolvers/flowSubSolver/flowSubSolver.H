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
    Foam::flowSubSolver

Description
    Base class for flow solvers

    \htmlinclude flowSubSolver_flowSubSolver.html

SourceFiles
    flowSubSolver.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2023

\*---------------------------------------------------------------------------*/

#ifndef flowSubSolver_H
#define flowSubSolver_H

#include "Time.H"
#include "IOdictionary.H"
#include "fvModels.H"
#include "fvConstraints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class flowSubSolver Declaration
\*---------------------------------------------------------------------------*/

class flowSubSolver
{
    // Private data     
                 
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        flowSubSolver(const flowSubSolver&);

        //- Disallow default bitwise assignment
        void operator=(const flowSubSolver&);
        
protected:

    // Protected data
        
        //- Reference to mesh
        const Time& runTime_;

        // Input/options dictionary
        const dictionary flowSubSolverDict_;

    // Protected Member Functions
    
public: 

    //- Runtime type information
    TypeName("none");  
    
    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        flowSubSolver,
        dictionary,
        (
            const Time& runTime,
            const dictionary& dict
        ),
        (runTime, dict)
    );
        
    // Constructors

        //- Construct from mesh,dictionary
        flowSubSolver
        (
            const Time& runTime,
            const dictionary& dict
        ); 

        //- Return a pointer to a new flowSubSolver
        static autoPtr<flowSubSolver> New
        (
            const Time& runTime,
            const dictionary& dict
        );

    //- Destructor
    virtual ~flowSubSolver();


    // Member Functions 
        
        //- Initialise time step before the solver loop
        virtual void init()
        {}

        //- Solve the subsolver equations
        virtual void correct()
        {}

        //- Finalise timestep after the solver loop
        virtual void finalise()
        {}

        //- Return true if converged
        virtual bool converged()
        {
            return true;
        }
        
        virtual scalar nextDeltaT() const
        {
            return great;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
