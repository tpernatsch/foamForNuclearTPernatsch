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
    Foam::isomolarBubble

Description
    Compute bubble/droplet diameter based on pressure and temperature
    changes assuming constant moles within bubble/droplet (i.e. no mass 
    exchange) via ideal gas law

Author
    Stefan Radman (stefan.radman@alumni.epfl.ch; sradman@protonmail.com)

SourceFiles
    isomolarBubbleFluidDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef isomolarBubble_H
#define isomolarBubble_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidDiameterModel.H"

namespace Foam
{
namespace fluidDiameterModels
{

/*---------------------------------------------------------------------------*\
            Class isomolarBubble Declaration
\*---------------------------------------------------------------------------*/

class isomolarBubble
:
    public fluidDiameterModel
{
protected:

    //- Reference values for diameter, pressure, temperature
    const scalar d0_;
    const scalar p0_;
    const scalar T0_;

    const volScalarField& p_;

    const volScalarField& T_;

public:

    //- Runtime type information
    TypeName("isomolarBubble");

    // Constructors

        isomolarBubble
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~isomolarBubble(){}

    // Member Functions

        //- Return fraction of structure interfacial contact in
        //  contact with the dispersed fluid
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidDiameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
