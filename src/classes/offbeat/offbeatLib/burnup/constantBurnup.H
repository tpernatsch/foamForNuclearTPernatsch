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
    Derived from burnup mother class, this class reads the burnup field from
    the initial condition file.

    \htmlinclude burnup_constantBurnup.html

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
    TypeName("fromLatestTime");

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
