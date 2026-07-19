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
    Foam::thermalSubSolver

Description
    When using the `constant` thermal solver class, the temperature field `T`
    remains unchanged from its initial configuration throughout the simulation.

    \note
    While the `constant` solver keeps the `T` field fixed by default, it is possible 
    to alter the `T` field during the simulation by leveraging advanced OpenFOAM 
    features like `#codeStream#`, `funkySetFields`, or specific function objects.

	\note
	At the moment, the `constant` class is one of two possible choices to
	couple OFFBEAT to an external solver that provides the `T` field directly at
	every time-step. At every coupled iteration, simply substitute the `T` file
	generated with the external solver in the starting time folder and (re)start the
	OFFBEAT simulation.

    \warning
    For backward-compatibility, `fromLatestTime` can be still used as the typename
    for the `constant` thermal solver. However, its use is deprecated.
	
________________________________________________________________________________

\par Options	
	As a `thermalSubSolver` class, `constant` **allows** the user to specify
	a few additional parameters in the `thermalOptions` sub-dictionary located 
	within the `solverDict`.
    
    Parameters in `thermalOptions`:
	
	- <b>`heatFluxSummary`</b> - It activates the calculation of the heat flowing out from
	every surface of the domain, printing this information on the terminal. **It is
	set to false by default**.
	- <b>`calculateEnthalpy`</b> - It activates the calculation of the enthalpy deposited
	in the fuel materials (useful for RIA cases). **It is set to false by default**.
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`constant` thermal class:
	
	\verbatim
        thermalSolver constant;

        thermalOptions
        {
            // Print summary of boundary-heat fluxes at each iteration
            heatFluxSummary false;

            // Activate calculation of enthalpy - false by default
            calculateEnthalpy false;
        }
	\endverbatim

SourceFiles
    thermalSubSolver.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef thermalSubSolver_H
#define thermalSubSolver_H

#include "materials.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class thermalSubSolver Declaration
\*---------------------------------------------------------------------------*/

class thermalSubSolver
{
    // Private data     
        
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        thermalSubSolver(const thermalSubSolver&);

        //- Disallow default bitwise assignment
        void operator=(const thermalSubSolver&);
        
protected:

    // Protected data
        
        //- Reference to mesh
        const fvMesh& mesh_;

        // Non const reference to the materials class
        materials& mat_;

        //- Input/options dictionary
        const dictionary thermalDict_;

        //- Temperature [K]
        volScalarField T_;

        //- Gradient of temperature field [K/m]
        volVectorField gradT_;
        
        //- Flag to determine if heat flux should be print to log and heatFlux 
        //- written to disk
        Switch heatFluxSummary_;

        //- Heat flux [W/m2]
        autoPtr<volScalarField> heatFlux_;
        autoPtr<volVectorField> heatFluxVec_;

        //- Switch to activate calculation of enthalpy
        Switch calculateEnthalpy_;

        //- Specific enthalpy [J/kg]
        autoPtr<volScalarField> enthalpy_;

        //- Residual of latest inner iteration
        scalar residual_;

        //- Residual of initial inner iteration 
        //- (first iteration of current inner iteration loop)
        scalar initialResidual_;  
        
        //- A different residual, based on the rate of change of T_
        scalar relResidual_;
        
        //- Absolute error, T-T.prevIter
        scalar absErr_;

    // Protected Member Functions
    
        //- Compute heat flux at the boundary
        void calcHeatFlux(const surfaceScalarField& phi);
    
        //- Compute heat flux at the boundary
        void calcEnthalpy();

        //- Compute temperature gradient 
        void computeGradT();

        //- Read and return number of inner correctors
        int nCorrectors(const word);
        
public: 

    //- Runtime type information
    TypeName("constant");  
    
    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        thermalSubSolver,
        dictionary,
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& thermalDict
        ),
        (mesh, mat, thermalDict)
    );
        
    // Constructors

        //- Construct from mesh, materials and dictionary
        thermalSubSolver
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& thermalDict
        ); 

        //- Return a pointer to a new thermalSolver
        static autoPtr<thermalSubSolver> New
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& solverDict
        );

    //- Destructor
    virtual ~thermalSubSolver();


    // Member Functions 

        //- Update temperature by solving the heat balance equation
        virtual void correct();

        //- Return true if converged given explicit residual values
        virtual bool converged
        (
            const scalar& residual,
            const scalar& relResidual,
            const scalar& absErr
        );

        //- Return true if converged using stored initial residuals
        virtual bool converged();

    //- Access functions     

        //- Return a const reference to the density through the material class
        inline const volScalarField& rho()
        {
            return mat_.rho();
        }

        //- Return a reference to the temperature field
        inline volScalarField& T()
        {
            return T_;
        }

        //- Return a const reference to the temperature field
        inline const volScalarField& T() const
        {
            return T_;
        }
        
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

        //- Return the maximum allowable timestep
        virtual scalar nextDeltaT() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
