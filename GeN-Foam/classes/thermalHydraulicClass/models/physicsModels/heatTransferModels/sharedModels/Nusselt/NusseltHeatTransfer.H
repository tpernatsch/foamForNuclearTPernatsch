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
    Foam::heatTransferModels::Nusselt

Description
    Heat transfer coefficient computed from Nusselt via a correlation in the
    form: A + B*Re^(C)*Pr^(D).

SourceFiles
    NusseltHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef Nusselt_H
#define Nusselt_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                           Class Nusselt Declaration
\*---------------------------------------------------------------------------*/

class Nusselt
:
    public heatTransferModel
{
private:

    //- Ref to bulk fluid
    const fluid& fluid_;

    //- Ref to Reynolds of pair
    const volScalarField& Re_;

    //- Ref to characteristic dimension of pair
    const volScalarField& Dh_;

    //- Coefficients in Nu = A_ + B_*Re^(C_)*Pr^(D_)
    const scalar A_;
    const scalar B_;
    const scalar C_;
    const scalar D_;

    //- For fluid-structure, marker field to set htc to 0 outside structure
    const scalarField& cellField_;

public:

    TypeName("Nusselt");


    //- Constructors

        //- Construct for a fluid-fluid pair
        Nusselt
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair,
            const word nameBulk,
            const word nameInterface
        );

        //- Construct for a fluid-structure pair
        Nusselt
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair,
            const word region
        );

        //- Destructor
        virtual ~Nusselt(){}


    //- Member Functions

        //- Drag coefficient
        virtual tmp<volScalarField> htc() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //