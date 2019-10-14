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
    Foam::fluidInterfacialAreaModels::Schor

Description
    Interfacial area model by Schor with hard-coded dependency for a 
    two-regimes map model. An interpolation regions is forcefully
    added between alpha1 and alpha2 to eliminate discontinuities.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    Schor.C

DictEntryExample ("" is the placeholder for the variable to be inserted)

    fluidInterfacialAreaModel
    {
        type        "Schor";
        alphaName   ""; // name of the vapour phase
        P           ""; // Pin pitch
        D           ""; // Pin diameter
    }

\*---------------------------------------------------------------------------*/

#ifndef Schor_H
#define Schor_H

#include "fluidInterfacialAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace fluidInterfacialAreaModels
{

/*---------------------------------------------------------------------------*\
                           Class Schor Declaration
\*---------------------------------------------------------------------------*/

class Schor
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

    //- P/D ratio
    const dimensionedScalar PD_;

    //- A constant
    const dimensionedScalar A_;

    //- PI
    const scalar PI_;

    //- Interpolation low alpha
    const dimensionedScalar alpha1_;

    //- Interpolation high alpha
    const dimensionedScalar alpha2_;

    //- Interpolation low iA
    const dimensionedScalar iA1_;

    //- Interpolation high iA
    const dimensionedScalar iA2_;


public:

    //- Runtime type information
    TypeName("Schor");


    // Constructors

        //- Construct for a fluid-fluid pair
        Schor
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid1,
            const fluid& fluid2
        );

    //- Destructor
    virtual ~Schor();


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