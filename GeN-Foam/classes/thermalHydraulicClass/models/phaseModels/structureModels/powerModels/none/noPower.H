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
    Foam::structureModels::none

Description

    Describes a structure with a constant (outer) temperature.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    none.C

\*---------------------------------------------------------------------------*/

#ifndef none_H
#define none_H

#include "powerModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerModels
{

/*---------------------------------------------------------------------------*\
                           Class none Declaration
\*---------------------------------------------------------------------------*/

class none
:
    public powerModel
{
protected :

public:

    //- Runtime type information
    TypeName("none");

    //- Constructors

        none
        (
            const structureModel& structure,
            const dictionary& dict,
            const word name
        );

        //- Destructor
        virtual ~none(){}

    //- Member Functions

        //- Return surface temperature (value does not matter, iA is set to 0,
        //  so no heat transfer will take place)
        virtual tmp<volScalarField> T() const;
        
        //- Do nothing
        virtual void correct
        (
            const volScalarField& HT,
            const volScalarField& H
        ) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //