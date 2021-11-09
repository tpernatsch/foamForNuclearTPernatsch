/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::saturationModel

Description

SourceFiles
    saturationModel.C
    newSaturationModel.C

\*---------------------------------------------------------------------------*/

#ifndef saturationModel_H
#define saturationModel_H

#include <map>
#include "fluid.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;
class phaseChangeModel;

/*---------------------------------------------------------------------------*\
                         Class saturationModel Declaration
\*---------------------------------------------------------------------------*/

class saturationModel
:
    public IOdictionary
{
protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to liquid
    const fluid& liquid_;

    //- Ref to vapour
    const fluid& vapour_;

    //- Pointer to a function (of the scalar value...(const label&) const)
    scalar (Foam::saturationModel::*funcPtr_)(const label&) const;

    //- Map that associates a function ptr to a word
    std::map<word, scalar (Foam::saturationModel::*)(const label&) const> map_;

public:

    TypeName("saturationModel");

    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationModel,
        saturationModels,
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        ),
        (pcm, dict, objReg)
    );

    //- Constructors

        saturationModel
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~saturationModel(){}

    //- Selectors

        static autoPtr<saturationModel> New
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

    //- Member Functions

        virtual scalar valueTSat(const label& celli) const = 0;

        virtual scalar valuePSat(const label& celli) const = 0;

        virtual scalar valueLnPSat(const label& celli) const = 0;

        virtual scalar valuePSatPrime(const label& celli) const = 0;

        //-
        virtual void correctField
        (
            volScalarField& field, 
            const word& func
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
