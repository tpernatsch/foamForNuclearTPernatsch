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
    Foam::porosityTransportSolver

Description
    Daugther class of transportSolver class that solves for the redistribution
    of porosity. 

    \htmlinclude transportSolver_porosityTransportSolver.html

SourceFiles
    porosityTransportSolver.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n

\contribution
    A. Scolaro, C. Fiorina - EPFL
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef porosityTransportSolver_H
#define porosityTransportSolver_H

#include "transportSolver.H"
#include "fuelMaterial.H"
#include "poreVelocityModel.H"
#include "layerAdditionRemovalPolyTopoChanger.H"
#include "displacementLaplacianFvMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class porosityTransportSolver Declaration
\*---------------------------------------------------------------------------*/

class porosityTransportSolver
:
public transportSolver
{
    // Private data
        
        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to temperature gradient field
        const volVectorField& gradT_;

        //- Constance reference to pore velocity
        const volVectorField& poreVelocity_;

        //- Ref to porosity field
        volScalarField& porosity_;

        //- Pore diffusion coefficient
        volScalarField d_;

        //- Field central void (=1 if void, =0 elsewhere)
        volScalarField centralVoid_;

        //- Field restructured region (=1 if restructured, =0 elsewhere)
        volScalarField restructuredRegion_;

        //- Threshold for central void field
        volScalarField asFabricatedP_;

        //- Dictionary for porosity options
        dictionary porosityDict_;

        //- Switch for diffusive component
        const bool diffusiveTerm_;

        //- Switch for advection term type
        const bool boundedAdvectionTerm_;

        //- Switch for mesh update
        const bool updateMesh_;

        //- Threshold for central void field
        scalar voidThreshold_;

        //- Criterion for time-stepping based on max Courant number
        scalar maxCoAllowed_;

    // Private Member Functions
        
        //- Set diffusion coefficient
        virtual void setDiffCoeff();

        //- Set centralVoid and restructuredRegion fields
        virtual void updateRestructuringFields();

        //- Disallow default bitwise copy construct
        porosityTransportSolver(const porosityTransportSolver&);

        //- Disallow default bitwise assignment
        void operator=(const porosityTransportSolver&);
     

public:

    //- Runtime type information
    TypeName("porosityTransport");

    // Constructors

        //- Construct from mesh, materials and dictionary
        porosityTransportSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        );
    
    //- Destructor
    ~porosityTransportSolver()
    {}

    // Access Functions

        // Return name of solving variable
        virtual const word& fieldName()
        {
            return porosity_.name();
        };

    // Member Functions
        
        //- Correct
        virtual void correct();

        //- Next time step according to solver criterion
        virtual scalar nextDeltaT();

        //- Update mesh 
        virtual void updateMesh()
        {
            // not implemented
        }
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
