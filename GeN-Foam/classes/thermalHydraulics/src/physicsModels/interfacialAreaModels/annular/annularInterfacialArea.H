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
    Interfacial area model that assumes that the dispersed phase flow is
    annular (i.e. cylindrical topology)

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    annular.C

\*---------------------------------------------------------------------------*/

#ifndef annularInterfacialArea_H
#define annularInterfacialArea_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfacialAreaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace interfacialAreaModels
{

/*---------------------------------------------------------------------------*\
            Class annular Declaration
\*---------------------------------------------------------------------------*/

class annular
:
    public interfacialAreaModel
{
protected:

    const volScalarField& dispersed_;
    const volScalarField& continuous_;
    const volScalarField& DhDispersed_;

    scalar cutoffAlpha_;

public:

    //- Runtime type information
    TypeName("annular");

    // Constructors

        annular
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~annular(){}

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
