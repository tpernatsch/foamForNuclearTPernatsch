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
    Foam::CHFModel::constantCHF

Description
    Gives a constant Critical Heat Flux provided in Phase Properties.

Author
    Gauthier Lazare (gauthier.lazare@gmail.com)

SourceFiles
    constantCHF.C

\*---------------------------------------------------------------------------*/

#ifndef constantCHF_H
#define constantCHF_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CHFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace criticalHeatFluxModels
{

/*---------------------------------------------------------------------------*\
            Class constantCHF Declaration
\*---------------------------------------------------------------------------*/

class constantCHF
:
    public CHFModel
{
protected:

    //- Constant Value of the CHF
    scalar value_;

public:

    //- Runtime type information
    TypeName("constantCHF");

    // Constructors

        constantCHF
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~constantCHF(){}

    // Member Functions

        //- Return CHF
        virtual scalar value
        (
            const label& celli
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace criticalHeatFluxModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
