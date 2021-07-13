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
    Foam::fluidDiameterModel

Description
    Class to handle fluid diameter

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    fluidDiameterModel.C
    newfluidDiameterModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidDiameterModel_H
#define fluidDiameterModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fluid.H"

namespace Foam
{

class FSPair;

/*---------------------------------------------------------------------------*\
            Class fluidDiameterModel Declaration
\*---------------------------------------------------------------------------*/

class fluidDiameterModel
:
    public IOdictionary
{
protected:

    const fvMesh& mesh_;

    const FSPair& pair_;

    const fluid& fluid_;

    const scalarField& dispersion_;

    //- Structure hydraulic diameter
    const volScalarField& Dhs_;

    //- Fluid diameter
    const volScalarField& Dh_;

public:

    //- Runtime type information
    TypeName("fluidDiameterModel");

    //- Declare runtime construction
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidDiameterModel,
            fluidDiameterModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    // Constructors

        fluidDiameterModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~fluidDiameterModel();

    // Selectors

        static autoPtr<fluidDiameterModel> New
        (
            const fluid& fluid,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        static autoPtr<fluidDiameterModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    // Member Functions

        //- Return diameter
        virtual scalar value(const label& celli) const = 0;

        //-
        void correctField(volScalarField& Dh) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
