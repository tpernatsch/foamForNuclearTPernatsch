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
    Foam::constantTabulatedAxialProfile

Description
    Axial profile defined using an constant interpolation table as a
    function of the axial position.
    Both "step" and "linear" interpolation methods are supported (see 
    Foam::InterpolateTable).
    The supplied profile must be normalised to 1.

SourceFiles
    constantTabulatedAxialProfile.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, E. Brunetto, C. Fiorina - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE 
    LAUSANNE, Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantTabulatedAxialProfile_H
#define constantTabulatedAxialProfile_H

#include "axialProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantTabulatedAxialProfile Declaration
\*---------------------------------------------------------------------------*/

class constantTabulatedAxialProfile
:
    public axialProfile
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantTabulatedAxialProfile(const constantTabulatedAxialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const constantTabulatedAxialProfile&);

public:

    //- Runtime type information
    TypeName("constantTabulated");

    
    // Constructors

        //- Construct from dictionary
        constantTabulatedAxialProfile
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            scalar zMin,
            scalar zMax 
        );

        
    //- Destructor
    virtual ~constantTabulatedAxialProfile();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
