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
    Foam::dragModels::Darcy

Description
    Drag coefficient in the form coeff*Re^exp. The coeff and exp are specified
    in vector form, one per local component (i.e. in the local structure frame
    of reference).

SourceFiles
    Darcy.C

\*---------------------------------------------------------------------------*/

#ifndef Darcy_H
#define Darcy_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class Darcy Declaration
\*---------------------------------------------------------------------------*/

class Darcy
:
    public dragModel
{
private:

    const FSPair& FSPair_;

    //-
    dimensionedVector coeffs_;

    //-
    dimensionedVector exps_;

public:

    TypeName("Darcy");

    //- Constructors

        //- Construct for a fluid-structure pair
        Darcy
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair,
            const word region
        );

        //- Destructor
        virtual ~Darcy(){}


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