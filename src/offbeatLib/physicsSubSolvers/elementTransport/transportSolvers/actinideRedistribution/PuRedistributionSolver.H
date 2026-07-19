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
    Foam::PuRedistributionSolver

Description
    Daugther class of transportSolver class that solves for the redistribution
    of Pu.
    
    \htmlinclude transportSolver_PuRedistributionSolver.html

SourceFiles
    PuRedistributionSolver.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n

\contribution
    A. Scolaro, C. Fiorina - EPFL
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef PuRedistributionSolver_H
#define PuRedistributionSolver_H

#include "transportSolver.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class PuRedistributionSolver Declaration
\*---------------------------------------------------------------------------*/

class PuRedistributionSolver
:
public transportSolver
{
    // Private data
        
        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to temperature gradient field
        const volVectorField& gradT_;

        //- Plutonium form factor
        volScalarField PuFormFactor_;

        //- Thermal diffusion coefficient Pu
        volScalarField diffPu_;

        //- Molar heat of transport [J/mol]
        const dimensionedScalar Q_;

        //- Universal gas constant [J/K/mol]
        const dimensionedScalar R_;

        //- Ref to porosity field
        const volScalarField& porosity_;

        //- Pores diameter [m]
        const dimensionedScalar poreDiameter_;

        //- Pores thickness  [m]
        const dimensionedScalar poreThickness_;

        //- Coefficient in pore migration term [1/K]
        const dimensionedScalar A_;

        //- Switch for pore migration term
        bool poreMigrationOn_;

    // Private Member Functions

        //- Initialize Pu concentration
        virtual void initializePuConcentration();

        //- Update diffusion coefficient and pore velocity
        virtual void updateCoefficients();

        //- Redistribute all Pu isotoped applying form factor
        virtual void redistributeIsotopes();

        //- Disallow default bitwise copy construct
        PuRedistributionSolver(const PuRedistributionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const PuRedistributionSolver&);
     

public:

    //- Runtime type information
    TypeName("PuRedistribution");

    // Constructors

        //- Construct from mesh, materials and dictionary
        PuRedistributionSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        );
    
    //- Destructor
    ~PuRedistributionSolver()
    {}

    // Access Functions

        // Return name of solving variable
        virtual const word& fieldName()
        {
            return PuFormFactor_.name();
        };

    // Member Functions
        
        //- Correct
        virtual void correct();

        //- Next time step according to solver criterion
        virtual scalar nextDeltaT()
        {
            //Not Implemented
            return GREAT;
        };

        virtual void updateMesh(){};
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
