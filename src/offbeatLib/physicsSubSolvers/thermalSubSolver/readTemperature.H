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
    Foam::readTemperature

Description
	The `readTemperature` thermal solver class in OFFBEAT reads the temperature `T`
	field from the time-step folders stored in the `constant/timeDependentInput/`
	sub-folder. If the `T` file is present in the corresponding time folder for a
	given time step, the temperature field will be set accordingly. Otherwise, the
	`T` field remains constant.
	
	\note
	The `readTemperature` class provides a second convenient way to (explicitly)
	couple OFFBEAT with an external solver that separately generates (before running
	OFFBEAT) the `T` field. To use this class, simply place the `T` files generated
	by the external solver in the corresponding time folders within the
	`constant/timeDependentInput/` sub-folder, and OFFBEAT will read and use that
	temperature field.

    \warning
    This class might be redundant, as the same functionalities might be achieved 
    using standard OpenFOAM functionObjects.
	
________________________________________________________________________________

\par Options	
	As a `thermalSubSolver` class, `readTemperature` **allows** the user to specify
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
	`readTemperature` thermal class:
	
	\verbatim
        thermalSolver readTemperature;

        thermalOptions
        {
            // Print summary of boundary-heat fluxes at each iteration
            heatFluxSummary false;

            // Activate calculation of enthalpy - false by default
            calculateEnthalpy false;
        }
	\endverbatim

SourceFiles
    readTemperature.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef readTemperature_H
#define readTemperature_H

#include "thermalSubSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class readTemperature Declaration
\*---------------------------------------------------------------------------*/

class readTemperature
:
public thermalSubSolver
{
    // Private data    

    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        readTemperature(const readTemperature&);

        //- Disallow default bitwise assignment
        void operator=(const readTemperature&);
     

public:

    //- Runtime type information
    TypeName("readTemperature");

    // Constructors

        //- Construct from mesh, materials and dictionary
        readTemperature
        (
            const fvMesh& mesh,
            materials& mat,
            const dictionary& thermalDict
        );
    
    //- Destructor
    ~readTemperature()
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
