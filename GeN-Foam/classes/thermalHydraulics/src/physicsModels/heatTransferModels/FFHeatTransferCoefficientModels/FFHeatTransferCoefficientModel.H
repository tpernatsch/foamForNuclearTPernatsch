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
    Foam::FFHeatTransferCoefficientModel

Description
    Class to manage the fluid-structure (i.e. fluid-wall) heat transfer
    coefficient

Author
    Stefan Radman (sradman@protonmail.com; stefan.radman@epfl.ch)

SourceFiles
    FFHeatTransferCoefficientModel.C
    newFFHeatTransferCoefficientModel.C

\*---------------------------------------------------------------------------*/

#ifndef FFHeatTransferCoefficientModel_H
#define FFHeatTransferCoefficientModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
                           Class FFHeatTransferCoefficientModel Declaration
\*---------------------------------------------------------------------------*/

class FFHeatTransferCoefficientModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FFPair& pair_;

    //- Ref to bulk fluid
    const fluid& bulkFluid_;

    //- Ref to other fluid
    const fluid& otherFluid_;

public:

    TypeName("FFHeatTransferCoefficientModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            FFHeatTransferCoefficientModel,
            FFHeatTransferCoefficientModels,
            (
                const FFPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        FFHeatTransferCoefficientModel
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~FFHeatTransferCoefficientModel(){}

    //- Selectors

        static autoPtr<FFHeatTransferCoefficientModel> New
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        virtual scalar value(const label& celli) const = 0;

        virtual void correctField(volScalarField& htc) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
