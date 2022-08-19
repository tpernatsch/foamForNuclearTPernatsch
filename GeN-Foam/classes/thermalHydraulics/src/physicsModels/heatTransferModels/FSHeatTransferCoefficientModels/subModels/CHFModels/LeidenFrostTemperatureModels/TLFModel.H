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
    Foam::TLFModel

Description
    Model for the prediction of the Leidenfrost temperature (i.e. T_min) for 
    which the heat flux is minimal. After this temperature, we are in post CHF Regime

Author
    Gauthier Lazare (gauthier.lazare@gmail.com)

SourceFiles
    TLFModel.C

\*---------------------------------------------------------------------------*/

#ifndef TLFModel_H
#define TLFModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structure.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

//-
class FSPair;

/*---------------------------------------------------------------------------*\
                           Class TLFModel Declaration
\*---------------------------------------------------------------------------*/

class TLFModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FSPair& pair_;

    //- Ref to Tsat (i.e. the interface T)
    const volScalarField& Tsat_;

    //- Ref to pressure 
    const volScalarField& p_;

public:

    TypeName("TLFModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            TLFModel,
            LeidenFrostTemperatureModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        TLFModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~TLFModel(){}

    //- Selectors

        static autoPtr<TLFModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions
        virtual scalar value
        (
            const label& celli
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
