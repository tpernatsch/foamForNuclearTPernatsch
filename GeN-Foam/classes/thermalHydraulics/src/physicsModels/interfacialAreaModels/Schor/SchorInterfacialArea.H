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

Description

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    Schor.C

\*---------------------------------------------------------------------------*/

#ifndef SchorInterfacialArea_H
#define SchorInterfacialArea_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfacialAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialAreaModels
{

/*---------------------------------------------------------------------------*\
            Class Schor Declaration
\*---------------------------------------------------------------------------*/

class Schor
:
    public interfacialAreaModel
{
protected:

    //- Phase fraction of vapour phase
    const fluid& vapour_;

    //- Pin diameter
    const scalar D_;

    //- Pin pitch
    const scalar P_;

    //- P/D ratio
    const scalar PD_;

    //- A constant
    const scalar A_;

    //- PI
    const scalar PI_;

    //- Interpolation low alpha
    const scalar alpha1_;

    //- Interpolation high alpha
    const scalar alpha2_;

    //-
    const scalar deltaAlpha21_;

    //- Cutoff alpha
    const scalar alpha3_;

    //- Interpolation low iA
    const scalar iA1_;

    //- Interpolation high iA
    const scalar iA2_;

    //- Min allowable iA at high alpha
    const scalar iA3_;

public:

    //- Runtime type information
    TypeName("Schor");

    // Constructors

        Schor
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~Schor(){}

    // Member Functions

        //- Return fraction of structure interfacial contact in
        //  contact with the dispersed fluid
        virtual scalar value(const label& celli) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace interfacialAreaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
