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
    Foam::isothermal

Description
    Compute bubble/droplet diameter based on pressure changes assuming
    constant temperature and mass within bubble/droplet (i.e. no mass exchange)
    via ideal gas law

Author
    Stefan Radman (stefan.radman@alumni.epfl.ch; sradman@protonmail.com)

SourceFiles
    isothermalBubbleFluidDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef isothermal_H
#define isothermal_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidDiameterModel.H"

namespace Foam
{
namespace fluidDiameterModels
{

/*---------------------------------------------------------------------------*\
            Class isothermal Declaration
\*---------------------------------------------------------------------------*/

class isothermal
:
    public fluidDiameterModel
{
protected:

    //- Reference values
    const scalar d0_;
    const scalar p0_;

    const volScalarField& p_;

public:

    //- Runtime type information
    TypeName("isothermal");

    // Constructors

        isothermal
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~isothermal(){}

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
