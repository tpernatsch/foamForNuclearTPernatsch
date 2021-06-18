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
    Foam::suppressionFactorModel

Description

Author
    Stefan Radman (sradman@protonmail.com; stefan.radman@epfl.ch)

SourceFiles
    suppressionFactorModel.C
    newSuppressionFactorModel.C

\*---------------------------------------------------------------------------*/

#ifndef suppressionFactorModel_H
#define suppressionFactorModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structure.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

//- Put them in a new FSHeatTransferCoefficientSubModels namespace? Eh, I don't 
//  know how handy this is

//-
class FSPair;

/*---------------------------------------------------------------------------*\
                           Class suppressionFactorModel Declaration
\*---------------------------------------------------------------------------*/

class suppressionFactorModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FSPair& pair_;

public:

    TypeName("suppressionFactorModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            suppressionFactorModel,
            suppressionFactorModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        suppressionFactorModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~suppressionFactorModel(){}

    //- Selectors

        static autoPtr<suppressionFactorModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        virtual scalar value(const label& celli) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
