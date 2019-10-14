/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::dragModels::NoKazimi

Description
    Imposes a constant drag coefficients between the phases.
    Howevere, the phase pair needs to be ordered. This is due to the fact
    that the momentum sink term itself (see in dragModel.C, function K())
    relies on continuous and disperesed data, which the present model 
    takes out of the way.
    By Stefan Radman (stefan.radman@epfl.ch)

SourceFiles
    NoKazimi.C

\*---------------------------------------------------------------------------*/

#ifndef NoKazimi_H
#define NoKazimi_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class phasePair;

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimi Declaration
\*---------------------------------------------------------------------------*/

class NoKazimi
:
    public dragModel
{

public:

    //- Runtime type information
    TypeName("NoKazimi");


    // Constructors

        //- Construct from a dictionary and a phase pair
        NoKazimi
        (
            const dictionary& dict,
            const phasePair& pair,
            const bool registerObject
        );


    //- Destructor
    virtual ~NoKazimi();


    // Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> Cd() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
