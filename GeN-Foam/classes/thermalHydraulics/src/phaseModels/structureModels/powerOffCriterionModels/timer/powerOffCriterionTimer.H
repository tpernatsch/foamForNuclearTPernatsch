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
    Foam::powerOffControlModels::timer

Description
    Turns power off if the field value maximum exceeds a provided value.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    timer.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    powerOffControlModel
    {
        type        "timer";
        value       "";   
    }


\*---------------------------------------------------------------------------*/

#ifndef timer_H
#define timer_H

#include "powerOffCriterionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace powerOffCriterionModels
{

/*---------------------------------------------------------------------------*\
                           Class timer Declaration
\*---------------------------------------------------------------------------*/

class timer
:
    public powerOffCriterionModel
{
private:

    const scalar t0_;

public:

    //- Runtime type information
    TypeName("timer");

    //- Constructors

        timer
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

    //- Destructor
    virtual ~timer();

    //- Member Functions

        //- Interfacial area
        virtual bool powerOffCriterion();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace powerOffControlModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //