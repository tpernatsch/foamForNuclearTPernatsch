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
    Foam::FFDragCoefficientModel

Description
    Class to manage the fluid-structure (i.e. fluid-wall) drag coefficient

Author
    Stefan Radman (sradman@protonmail.com; stefan.radman@epfl.ch)

SourceFiles
    FFDragCoefficientModel.C
    newFFDragCoefficientModel.C

\*---------------------------------------------------------------------------*/

#ifndef FFDragCoefficientModel_H
#define FFDragCoefficientModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
                           Class FFDragCoefficientModel Declaration
\*---------------------------------------------------------------------------*/

class FFDragCoefficientModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FFPair& pair_;

    //-
    const volScalarField& Re_;

public:

    TypeName("FFDragCoefficientModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            FFDragCoefficientModel,
            FFDragCoefficientModels,
            (
                const FFPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        FFDragCoefficientModel
        (
            const FFPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~FFDragCoefficientModel(){}

    //- Selectors

        static autoPtr<FFDragCoefficientModel> New
        (
            const FFPair& pair,
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
