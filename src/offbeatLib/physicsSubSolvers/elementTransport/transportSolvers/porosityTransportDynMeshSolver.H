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
    Foam::porosityTransportDynMeshSolver

Description
    Daugther class of transportSolver class that solves for the redistribution
    of porosity. 

    \htmlinclude transportSolver_porosityTransportDynMeshSolver.html

SourceFiles
    porosityTransportDynMeshSolver.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n

\contribution
    A. Scolaro, C. Fiorina - EPFL
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef porosityTransportDynMeshSolver_H
#define porosityTransportDynMeshSolver_H

#include "transportSolver.H"
#include "fuelMaterial.H"
#include "poreVelocityModel.H"
#include "layerAdditionRemovalPolyTopoChanger.H"
#include "displacementLaplacianFvMotionSolver.H"
#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class porosityTransportDynMeshSolver Declaration
\*---------------------------------------------------------------------------*/

class porosityTransportDynMeshSolver
:
public transportSolver
{
    // Private data

        //- Reference to dynamicFvMesh (casted from fvMesh)
        dynamicFvMesh& dynMesh_;

        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to temperature gradient field
        const volVectorField& gradT_;

        //- Constance reference to pore velocity
        volVectorField& poreVelocity_;

        //- Ref to porosity field
        volScalarField& porosity_;

        //- Pore diffusion coefficient
        volScalarField d_;

        //- Dictionary for porosity options
        dictionary porosityDict_;

        //- Switch for diffusive component
        const bool diffusiveTerm_;

        //- Switch for advection term type
        const bool boundedAdvectionTerm_;

        //- Switch for mesh update
        const bool updateMesh_;

        //- Inner hole patch IDs (useful for central void opening)
        labelList holePatchIDs_;

        //- Criterion for time-stepping based on max Courant number
        scalar maxCoAllowed_;

        //- Coalescence velocity parameters
        scalarList vCoalParameters_;

        //- Settings related to initialization of central hole position
        Switch moveHole_;
        vector holePosition_;
        bool holeMoved_;
        bool holeMoving_;
        
        //- Enable correction for the increasing inner radius in lower scale model
        Switch useRiCorrection_;

        //- Use direction pore velocity (inverse) for the expansion of central hole
        Switch useDirectionPV_;
        
    // Private Member Functions
        
        //- Set diffusion coefficient
        virtual void setDiffCoeff();

        //- Read and initialize hole patch list
        virtual void readHolePatches();

        //- Disallow default bitwise copy construct
        porosityTransportDynMeshSolver(const porosityTransportDynMeshSolver&);

        //- Disallow default bitwise assignment
        void operator=(const porosityTransportDynMeshSolver&);
     

public:

    //- Runtime type information
    TypeName("porosityTransportDynMesh");

    // Constructors

        //- Construct from mesh, materials and dictionary
        porosityTransportDynMeshSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        );
    
    //- Destructor
    ~porosityTransportDynMeshSolver()
    {}

    // Access Functions

        // Return name of solving variable
        virtual const word& fieldName()
        {
            return porosity_.name();
        };

    // Member Functions

        //- Compute current hole center
        virtual vector computeHoleCenter();

        //- Move position of initialized central hole
        virtual void moveInitialHole();

        //- Compute coalescence velocity for given patch ID
        virtual const scalarField coalescenceVelocity
        (
            const label patchID
        );

        //- Update mesh 
        virtual void updateMesh();
        
        //- Correct
        virtual void correct();

        //- Next time step according to solver criterion
        virtual scalar nextDeltaT();
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
