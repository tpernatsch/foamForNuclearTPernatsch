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
    Foam::fluidDiameterModels::isothermal

Description
    isothermal fluid diameter model. Adapted from twoPhaseEulerFoam

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    isothermal.C

DictEntryExample

    fluidDiameterModel
    {
        type        isothermal;
        d0          reference diameter for isothermal expansion;
        p0          reference pressure for isothermal expansion;
    }


\*---------------------------------------------------------------------------*/

#ifndef isothermal_H
#define isothermal_H

#include "fluidDiameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

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
private:

    //- Reference to system pressure
    const volScalarField p_;

    //- Reference diameter for the isothermal expansion
    dimensionedScalar d0_;

    //- Reference pressure for the isothermal expansion
    dimensionedScalar p0_;

public:

    //- Runtime type information
    TypeName("isothermal");

    //- Constructors

        isothermal
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid
        );

    //- Destructor
    virtual ~isothermal();

    //- Member Functions

        //- Interfacial area
        virtual tmp<volScalarField> Dh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidDiameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //