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
    Foam::solidConductionSolver

Description
    Class handling the solid heat conduction in presence of a heat source.

    \htmlinclude thermalSubSolver_solidConductionSolver.html

SourceFiles
    solidConductionSolver.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef solidConductionSolver_H
#define solidConductionSolver_H

#include "thermalSubSolver.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class solidConductionSolver Declaration
\*---------------------------------------------------------------------------*/

class solidConductionSolver
:
public thermalSubSolver
{
    // Private data   
        
        //- Reference to heatSource field
        const volScalarField* Q_;   

        //- Gas gap heat conductance
        volScalarField hGap_;

        // Private Member Functions

        //- Disallow default bitwise copy construct
        solidConductionSolver(const solidConductionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const solidConductionSolver&);
     

public:

    //- Runtime type information
    TypeName("solidConduction");

    // Constructors

        //- Construct from mesh, materials and dictionary
        solidConductionSolver
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& thermalOptDict
        );
    
    //- Destructor
    ~solidConductionSolver()
    {}

    // Member Functions
        
        //- Update temperature by solving the solid heat conduction equation
        virtual void correct();  
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
