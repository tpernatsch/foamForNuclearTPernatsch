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
    Foam::pimpleFlowSubSolver

Description
    Flow solver based using PIMPLE algorithm

    \htmlinclude flowSubSolver_pimpleFlowSubSolver.html

SourceFiles
    pimpleFlowSubSolver.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2023

\*---------------------------------------------------------------------------*/

#ifndef pimpleFlowSubSolver_H
#define pimpleFlowSubSolver_H

#include "flowSubSolver.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "singlePhaseTransportModel.H"
#include "momentumTransportModel.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "IOMRFZoneList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class pimpleFlowSubSolver Declaration
\*---------------------------------------------------------------------------*/

class pimpleFlowSubSolver
:
public flowSubSolver
{
    // Private data
    
        //- Pointer to mesh (if needed)
        autoPtr<fvMesh> meshPtr_;
        
        //- Reference to mesh
        fvMesh& mesh_;
        
        //- Pressure (Pa)
        volScalarField p_;
        
        //- Velocity (m/s)
        volVectorField U_;
        
        //- Surface mass fluxes
        surfaceScalarField phi_;
        
        //- Transport model
        singlePhaseTransportModel laminarTransport_;
        
        //-Turbulence model
        autoPtr<incompressible::momentumTransportModel> turbulence_;

        //- PIMPLE solver controls
        pimpleControl pimple_;
        
        //- Reference pressure
        pressureReference pressureReference_;

        //- Continuity Errors
        scalar cumulativeContErr_;

        //- Pointer to the demand driven fvModels MeshObject
        mutable Foam::fvModels fvModels_;

        //- Pointer to the demand driven fvConstraints MeshObject
        mutable Foam::fvConstraints fvConstraints_;
        
    // Protected Member Functions
        
        //- Solve the momentum equation
        tmp<fvVectorMatrix> solveMomentum();

        //- Solve the pressure correction equation
        void correctPressure(fvVectorMatrix& UEqn);
        
        //- Calculate and print the continuity errors
        void continuityErrors();

    // Private Member Functions
        
        //- Helper function to return a reference to the mesh or construct as needed
        fvMesh& getMesh(const Time& runTime);

        //- Disallow default bitwise copy construct
        pimpleFlowSubSolver(const pimpleFlowSubSolver&);

        //- Disallow default bitwise assignment
        void operator=(const pimpleFlowSubSolver&);
        

public:

    //- Runtime type information
    TypeName("pimple");

    // Constructors

        //- Construct from mesh, materials and dictionary
        pimpleFlowSubSolver
        (
            const Time& runTime,
            const dictionary& dict
        );
    
    //- Destructor
    ~pimpleFlowSubSolver()
    {}

    // Member Functions

        //- Return the fvModels that are created on demand
        inline Foam::fvModels& fvModels() const
        {
            return fvModels_;
        }
        
        //- Return the fvConstraints that are created on demand
        inline Foam::fvConstraints& fvConstraints() const
        {
            return fvConstraints_;
        }

        //- Updates before the PIMPLE loop
        virtual void init();

        //- Updates within the PIMPLE loop
        virtual void correct();

        //- Updates after the PIMPLE loop
        virtual void finalise();

        //- Return true if converged
        virtual bool converged();
        
        virtual scalar nextDeltaT() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
