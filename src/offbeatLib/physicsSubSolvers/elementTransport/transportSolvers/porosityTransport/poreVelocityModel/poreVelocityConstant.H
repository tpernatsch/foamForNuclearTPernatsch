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
    Foam::poreVelocityConstant

Description
    By selecting this constant model the velocity is read from the previous time
    step folder. Useful for user-defined velocity functions (e.g. by means of 
    codedStream-type initial conditions).

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UPuO2;

            ...

            poreVelocityModel     constant;
        }

        ...
    }
    \endverbatim

SourceFiles
    poreVelocityConstant.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef poreVelocityConstant_H
#define poreVelocityConstant_H

#include "poreVelocityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poreVelocityConstant Declaration
\*---------------------------------------------------------------------------*/

class poreVelocityConstant
:
    public poreVelocityModel
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        poreVelocityConstant(const poreVelocityConstant&);

        //- Disallow default bitwise assignment
        void operator=(const poreVelocityConstant&);

protected:
 

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from mesh, materials, and dictionary
        poreVelocityConstant
        (
            const fvMesh& mesh, 
            const dictionary& dict,
            const word defaultModel
        );


    //- Destructor
    ~poreVelocityConstant();

    // Access functions    

    // Member Functions

        //- Correct velocity 
        virtual void correct( const labelList& addr );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
