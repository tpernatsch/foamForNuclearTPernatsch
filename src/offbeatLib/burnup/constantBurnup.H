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
    Foam::constantBurnup

Description
	The `constant` burnup class in OFFBEAT allows you to set the burnup field
	in a file located in the starting time folder. If the `Bu` file is present in the
	starting time folder, the internal field and boundary conditions are set
	accordingly. If the `Bu`file is not present, the burnup field is set by degault to 
    0 MWd/MT\f$_{oxide}\f$ in each cell, with zeroGradient BCs in all non-empty/-wedge
	patches. oNte that.
    
    When selecting this burnup model, the burnup field will not evolve over time and
    will remain as defined in the starting time folder (or as set by default by OFFBEAT).
	
	\note
	At the moment, the `constant` model is the correct burnup model to use if
	you are coupling OFFBEAT to an external solver that provides the `Bu` field
	directly. At every coupled iteration, simply substitute the `Bu` file generated
	with the external solver in the starting time folder and (re)start the OFFBEAT
	simulation.

    \warning
    For backward-compatibility, `fromLatestTime` can be still used as the typename
    for the `constant` burnup model. However, its use is deprecated.
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`constant` burnup model
	
	\verbatim
	    burnup constant;
	\endverbatim

SourceFiles
    constantBurnup.C 

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantBurnup_H
#define constantBurnup_H

#include "burnup.H"
#include "volFields.H"
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantBurnup Declaration
\*---------------------------------------------------------------------------*/

class constantBurnup
:
    public burnup
{
    // Private data

protected:
    
    // Protected Member Functions
    
        //- Burnup in MWd/MT
        volScalarField Bu_;

        //- Maximum burnup increase for adjustable time step
        scalar maxBuIncrease_;
        
        //- Name used for the heatSource field
        word heatSourceName_;  
    
    // Protected Member Functions

        //- Disallow default bitwise copy construct
        constantBurnup(const constantBurnup&);

        //- Disallow default bitwise assignment
        void operator=(const constantBurnup&);

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from mesh, materials and dict
        constantBurnup
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& burnupDict
        );

        //- Destructor
        ~constantBurnup();


    // Member Functions

        //- Return a constant reference to the burnup field
        virtual const volScalarField& Bu() const
        {
            return Bu_;
        };
        
        //- Return next time increment according to burnup model
        virtual scalar nextDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
