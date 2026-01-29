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
    Foam::burnupFromPower

Description
	The `fromPower` burnup class can be used to track the evolution of the local
	burnup in MWd/MT\f$_{oxide}\f$.
	
	#### Formulation

	At each time step, the updated burnup is calculated as:

	$$
	\begin{aligned}
	Bu = Bu^0 + \frac{Q \cdot \Delta t}{1000 \cdot \rho}
	\end{aligned}
	$$
	
	where:
	
	- \f$Bu^0\f$ is the old local burnup value, i.e. at the end of the previous time
	step, in MWd/MT\f$_{oxide}\f$
	- \f$Q\f$ is the power density in W/m\f$^3\f$
	- \f$\rho\f$ is density in kg/m\f$^3\f$
	- \f$\Delta t\f$ is time step in seconds.
	
________________________________________________________________________________

\par Options
	The `fromPower` burnup model has at the moment only a single additional
	parameter that the user **might** specify in the `burnupOptions` sub-dictionary.

    Parameters in `burnupOptions`:

    - <b>`heatSourceName`</b> - The name of the heat source field. **By default this is
    set as `Q` in OFFBEAT**. This keyword could be useful if the power density field 
    (for instance calculted via a separate code) is provided via a different/external 
    heat source model with a different name.

________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`fromPower` burnup model
	
	\verbatim
        burnup fromPower;

        burnupOptions
        {
            // Name of the heat source field (set to "Q" by default)
            heatSourceName  "Q";
        }
	\endverbatim

SourceFiles
    burnupFromPower.C 

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef burnupFromPower_H
#define burnupFromPower_H

#include "constantBurnup.H"
#include "volFields.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class burnupFromPower Declaration
\*---------------------------------------------------------------------------*/

class burnupFromPower
:
    public constantBurnup
{
    // Private data

protected:
    
    //- Name used for the heatSource field
    word heatSourceName_;  
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        burnupFromPower(const burnupFromPower&);

        //- Disallow default bitwise assignment
        void operator=(const burnupFromPower&);

public:

    //- Runtime type information
    TypeName("fromPower");

    // Constructors

        //- Construct from mesh, materials and dict
        burnupFromPower
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& burnupDict
        );

        //- Destructor
        ~burnupFromPower();


    // Member Functions
    
        //- Update the burnup
        virtual void correct();

        //- Return error for this class
        virtual const scalarField SigmaA() const
        {
            FatalErrorIn("burnupFromPower::SigmaA()")
            << "The field SigmaA does not exist if a burnup depletion model is not "
            << "enabled. Either use the 'burnupLassman' burnup model or switch "
            << "the neutronics subsolver to 'fromLatestTime'."
            << exit(FatalError);

            return volScalarField::null(); // Dummy return field (avoid warning)
        };

        //- Return error for this class
        virtual const scalarField neutronD() const
        {
            FatalErrorIn("burnupFromPower::neutronD()")
            << "The field neutronD does not exist if a burnup depletion model is not "
            << "enabled. Either use the 'burnupLassman' burnup model or switch "
            << "the neutronics subsolver to 'fromLatestTime'."
            << exit(FatalError);

            return volScalarField::null(); // Dummy return field (avoid warning)
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
