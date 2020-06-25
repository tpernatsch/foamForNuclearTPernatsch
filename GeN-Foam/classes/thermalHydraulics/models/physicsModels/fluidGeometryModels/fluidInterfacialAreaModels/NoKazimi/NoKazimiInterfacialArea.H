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
    Foam::fluidInterfacialAreaModels::NoKazimiInterfacialArea

Description
    Interfacial area model by NoKazimiInterfacialArea with hard-coded 
    dependency for a  two-regimes map model. An interpolation regions is 
    added between alpha1 and alpha2 to eliminate discontinuities.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    NoKazimiInterfacialArea.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    fluidInterfacialAreaModel
    {
        type        "NoKazimiInterfacialArea";
        P           ""; // Pin pitch
        D           ""; // Pin diameter
    }

\*---------------------------------------------------------------------------*/

#ifndef NoKazimiInterfacialArea_H
#define NoKazimiInterfacialArea_H

#include "fluidInterfacialAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidInterfacialAreaModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimiInterfacialArea Declaration
\*---------------------------------------------------------------------------*/

class NoKazimiInterfacialArea
:
    public fluidInterfacialAreaModel
{
private:

    //- Phase fraction of vapour phase
    const volScalarField& alpha_;

    //- Pin diameter
    const dimensionedScalar D_;

    //- Pin pitch
    const dimensionedScalar P_;

    //- P/D
    const dimensionedScalar PD_;

    //- A constant
    const dimensionedScalar A_;

public:

    //- Runtime type information
    TypeName("NoKazimi");


    // Constructors

        //- Construct for a fluid-fluid pair
        NoKazimiInterfacialArea
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous
        );

    //- Destructor
    virtual ~NoKazimiInterfacialArea();


    // Member Functions

        //- Interfacial area
        virtual tmp<volScalarField> iA() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidInterfacialAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //