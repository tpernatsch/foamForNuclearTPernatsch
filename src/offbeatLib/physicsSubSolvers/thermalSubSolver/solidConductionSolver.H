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
	In the `solidConduction` solver, the temperature distribution is obtained from
	the solution of the heat diffusion or energy conservation equation.	

	#### Formulation
    
	For an arbitrary body with volume V, bounded by a surface S with outward normal
	n, this equation can be expressed as:
	
	$$
	\begin{aligned}
	\int \rho c_p \frac{dT}{dt} dt = \int_S \vec{q''} \cdot \vec{n} dS + \int_V Q dV
	\end{aligned}
	$$
	
	where:
	
	- \f$\rho\f$ is density kg/m\f$^3\f$
	- \f$c_p\f$ is the specific heat capacity in J/kg/K/
	- \f$Q\f$ is the volumetric heat source in W/m\f$^3\f$
	- \f$q''\f$ is the surface heat flux W/m\f$^2\f$
	
	Using Fourier's law, one can express the heat flux in terms of the temperature
	gradient as:
	
	$$
	\begin{aligned}
	\vec{q''} = k \nabla T
	\end{aligned}
	$$
	
	where \f$k\f$ is thermal conductivity in W/m/k.
	
	Thus, the final form of the equation (applying the Gauss-Green theorem) is:
	
	$$
	\begin{aligned}
	\int \rho c_p \frac{dT}{dt} dt = \int_V \nabla(k\nabla T) dV + \int_V Q dV
	\end{aligned}
	$$

________________________________________________________________________________

\par Options	
	As a `thermalSubSolver` class, `solidConduction` requires the user to specify a
	few additional parameters in the `thermalOptions` sub-dictionary located 
	within the `solverDict`.
    
    Parameters in `thermalOptions`:
	
	- <b>`heatFluxSummary`</b> - It activates the calculation of the heat flowing out from
	every surface of the domain, printing this information on the terminal. **It is
	set to false by default**.
	- <b>`calculateEnthalpy`</b> - It activates the calculation of the enthalpy deposited
	in the fuel materials (useful for RIA cases). **It is set to false by default**.
	
________________________________________________________________________________

Usage
	Here is an example of the `solverDict` to be used for activating the
	`solidConduction` thermal solver:
	
	\verbatim
		thermalSolver solidConduction;

		thermalOptions
		{
			// Print summary of boundary-heat fluxes at each iteration
			heatFluxSummary false;

			// Activate calculation of enthalpy - false by default
			calculateEnthalpy false;
		}
	\endverbatim

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
            const dictionary& thermalDict
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
