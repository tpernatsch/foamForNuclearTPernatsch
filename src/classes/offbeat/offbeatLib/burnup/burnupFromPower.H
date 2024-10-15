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
    Derived from constantBurnup, this class calculates the evolution of burnup 
    over time from the local power density.

    \htmlinclude burnup_burnupFromPower.html

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
