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
    Foam::dragModels::SchillerNaumann

Description
    SchillerNaumann correlation for the fluid-fluid drag coefficient. One among
    the most widely employed in multi-phase flows

Author
    Stefan Radman (stefan.radman@epfl.ch)

SourceFiles
    SchillerNaumann.C

\*---------------------------------------------------------------------------*/

#ifndef SchillerNaumann_H
#define SchillerNaumann_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class SchillerNaumann Declaration
\*---------------------------------------------------------------------------*/

class SchillerNaumann
:
    public dragModel
{
protected:

public:

    TypeName("SchillerNaumann");

    //- Constructors

        //- Construct from a dictionary and a phase pair
        SchillerNaumann
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair
        );

        //- Destructor
        virtual ~SchillerNaumann(){}

    //- Member Functions

        //- Drag coefficient
        virtual void correctKd(volTensorField& Kd) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
