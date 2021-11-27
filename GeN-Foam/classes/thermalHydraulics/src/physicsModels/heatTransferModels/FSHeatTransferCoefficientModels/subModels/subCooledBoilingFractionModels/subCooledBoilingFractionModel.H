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
    Foam::subCooledBoilingFractionModel

Description
    Model for the prediction of the sub-cooled bulk liquid temperature at 
    which bubble detachment occurs (only for sub-cooled boiling)

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    subCooledBoilingFractionModel.C
    newsubCooledBoilingFractionModel.C

\*---------------------------------------------------------------------------*/

#ifndef subCooledBoilingFractionModel_H
#define subCooledBoilingFractionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structure.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class FSPair;
class FFPair;

/*---------------------------------------------------------------------------*\
                           Class subCooledBoilingFractionModel Declaration
\*---------------------------------------------------------------------------*/

class subCooledBoilingFractionModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FSPair& pair_;

    //- Ptr to FFPair (fUcK eNcApSuLaTi1o0N!11!1)
    mutable const FFPair* FFPairPtr_;

public:

    TypeName("subCooledBoilingFractionModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            subCooledBoilingFractionModel,
            subCooledBoilingFractionModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        subCooledBoilingFractionModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~subCooledBoilingFractionModel(){}

    //- Selectors

        static autoPtr<subCooledBoilingFractionModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        //-
        void setPtr() const;

        //- Return value of the fraction of the pool boiling heat flux
        //  that is resulting in a net vapour generation
        virtual scalar value
        (
            const label& celli,
            const scalar& qi
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
