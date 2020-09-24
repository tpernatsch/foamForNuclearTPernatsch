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
    Foam::heatTransferModels::NoKazimiFFHeatTransfer

Description
    Heat transfer coefficient computed from Nusselt via a correlation in the
    form: A + B*Reynolds^(C)*Prandtl^(D).

Author
    Stefan Radman (stefan.radman@epfl.ch, stefanradman92@gmail.com)

SourceFiles
    NoKazimiFFHeatTransfer.C

\*---------------------------------------------------------------------------*/

#ifndef NoKazimiFFHeatTransfer_H
#define NoKazimiFFHeatTransfer_H

#include "heatTransferModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace heatTransferModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimiFFHeatTransfer Declaration
\*---------------------------------------------------------------------------*/

class NoKazimiFFHeatTransfer
:
    public heatTransferModel
{
private:

    //- Ref to bulk fluid
    const fluid& bulkFluid_;

    //- Ref to other fluid across the interface
    const fluid& otherFluid_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Ref to fluid-fluid interfacial area
    const volScalarField& iA_;

    //- Ref to fluid-fluid interfacial temperature
    const volScalarField& iT_;

    //- A pre-computed coefficients
    dimensionedScalar coeff_;

public:

    TypeName("NoKazimi");


    //- Constructors

        //- Construct for a fluid-fluid pair
        NoKazimiFFHeatTransfer
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair,
            const word nameBulk,
            const word nameInterface
        );

        //- Destructor
        virtual ~NoKazimiFFHeatTransfer(){}


    //- Member Functions

        //- Drag coefficient
        virtual void correctHtc(volScalarField& htc) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //