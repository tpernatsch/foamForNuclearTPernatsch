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
    Foam::neutronicsSubSolver

Description
    Parent class for neutronics subsolver. 

    \htmlinclude neutronicsSubSolver_neutronicsSubSolver.html

SourceFiles
    neutronicsSubSolver.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n

\contribution
    A. Scolaro, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef neutronicsSubSolver_H
#define neutronicsSubSolver_H

#include "burnup.H"
#include "fvMesh.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class neutronicsSubSolver Declaration
\*---------------------------------------------------------------------------*/

class neutronicsSubSolver
{
    // Private data     
        
    // Private Member Functions
        
    //- Disallow default bitwise copy construct
    neutronicsSubSolver(const neutronicsSubSolver&);

    //- Disallow default bitwise assignment
    void operator=(const neutronicsSubSolver&);
        
protected:
        
    //- Reference to mesh
    const fvMesh& mesh_;

    // Const reference to the burnup class
    const burnup& burnup_;

    // Reference to the neutronicsSubSolver option dictionary
    const dictionary& neutronicsOptDict_;

    //- Residual of latest inner iteration
    scalar residual_;
    
    //- Residual of initial inner iteration 
    scalar initialResidual_;  
    
    //- Relative residual
    scalar relResidual_;
    
    // Protected Member Functions
    
        
public: 

    //- Runtime type information
    TypeName("fromLatestTime");  
    
    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        neutronicsSubSolver,
        dictionary,
        (
            const fvMesh& mesh,
            const burnup& bu,
            const dictionary& neutronicsOptDict
        ),
        (mesh, bu, neutronicsOptDict)
    );
        
    // Constructors

    //- Construct from mesh, burnup and dictionary
    neutronicsSubSolver
    (
        const fvMesh& mesh,
        const burnup& bu,
        const dictionary& neutronicsOptDict
    ); 

    //- Return a pointer to a new neutronicsSubSolver
    static autoPtr<neutronicsSubSolver> New
    (
        const fvMesh& mesh,
        const burnup& bu,
        const dictionary& solverDict
    );

    //- Destructor
    virtual ~neutronicsSubSolver();


    // Member Functions 

    //- Update neutron flux 
    virtual void correct(){};

    //- Return true if converged
    virtual bool converged()
    {
        return true;
    }

    //- Access functions
    
    //- Return a const reference to residual_
    inline const scalar& residual() 
    {
        return residual_;
    };
    
    //- Return a const reference to initialResidual_
    inline const scalar& initialResidual() 
    {
        return initialResidual_;
    }; 
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
