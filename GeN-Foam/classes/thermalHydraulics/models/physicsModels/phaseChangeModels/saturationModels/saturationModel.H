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

#include "IOdictionary.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;

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

    //- Ref to fluid1
    const fluid& fluid1_;

    //- Ref to fluid2
    const fluid& fluid2_;

    //- Dimensioned ones needed for dimless conversions in log and exp funcs
    dimensionedScalar onePa_;
    dimensionedScalar oneK_;

public:

    TypeName("saturationModel");

    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        saturationModel,
        saturationModels,
        (
            const dictionary& dict, 
            const fluid& fluid1,
            const fluid& fluid2
        ),
        (dict, fluid1, fluid2)
    );

    //- Constructors

        saturationModel
        (
            const dictionary& dict, 
            const fluid& fluid1,
            const fluid& fluid2
        );

        //- Destructor
        virtual ~saturationModel(){}

    //- Selectors

        static autoPtr<saturationModel> New
        (
            const dictionary& dict, 
            const fluid& fluid1,
            const fluid& fluid2
        );

    //- Member Functions

        //- Saturation temperature at given pressure
        virtual tmp<volScalarField> Tsat
        (
            const volScalarField& p
        ) const = 0;

        //- Saturation pressure at given temperature
        virtual tmp<volScalarField> pSat
        (
            const volScalarField& T
        ) const = 0;

        //- Natural log of the saturation pressure at given temperature
        virtual tmp<volScalarField> lnPSat
        (
            const volScalarField& T
        ) const = 0;   

        //- dPressureSat/dTemperature at given temperature
        virtual tmp<volScalarField> pSatPrime
        (
            const volScalarField& T
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
