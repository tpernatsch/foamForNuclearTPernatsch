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
    Foam::AmRedistributionSolver

Description
    Daugther class of transportSolver class that solves for the redistribution
    of Am. 

    \htmlinclude transportSolver_AmRedistributionSolver.html
    
SourceFiles
    AmRedistributionSolver.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n

\contribution
    A. Scolaro, C. Fiorina - EPFL
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef AmRedistributionSolver_H
#define AmRedistributionSolver_H

#include "transportSolver.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class AmRedistributionSolver Declaration
\*---------------------------------------------------------------------------*/

class AmRedistributionSolver
:
public transportSolver
{
    // Private data
        
        //- Reference to temperature field
        const volScalarField& T_;

        //- Reference to temperature gradient field
        const volVectorField& gradT_;

        //- Plutonium form factor
        volScalarField AmFormFactor_;

        //- Thermal diffusion coefficient Am
        volScalarField diffAm_;

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

        //- Switch for pore migration
        bool poreMigrationOn_;

    // Private Member Functions

        //- Initialize Am concentration
        virtual void initializeAmConcentration();

        //- Update diffusion coefficient and pore velocity
        virtual void updateCoefficients();

        //- Disallow default bitwise copy construct
        AmRedistributionSolver(const AmRedistributionSolver&);

        //- Disallow default bitwise assignment
        void operator=(const AmRedistributionSolver&);
     

public:

    //- Runtime type information
    TypeName("AmRedistribution");

    // Constructors

        //- Construct from mesh, materials and dictionary
        AmRedistributionSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        );
    
    //- Destructor
    ~AmRedistributionSolver()
    {}

    // Access Functions

        // Return name of solving variable
        virtual const word& fieldName()
        {
            return AmFormFactor_.name();
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
