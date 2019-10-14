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
    Foam::dragModels::constant

Description
    Constant drag coefficient model

SourceFiles
    constant.C

\*---------------------------------------------------------------------------*/

#ifndef constant_H
#define constant_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class constant Declaration
\*---------------------------------------------------------------------------*/

class constant
:
    public dragModel
{
private:

    volTensorField Kd_;

public:

    TypeName("constant");

    //- Constructors

        //- Construct for a fluid-fluid pair
        constant
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair
        );

        //- Construct for a fluid-structure pair
        constant
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair,
            const word region
        );

        //- Destructor
        virtual ~constant(){}


    //- Member Functions

        //- Drag coefficient
        virtual tmp<volTensorField> Kd() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //