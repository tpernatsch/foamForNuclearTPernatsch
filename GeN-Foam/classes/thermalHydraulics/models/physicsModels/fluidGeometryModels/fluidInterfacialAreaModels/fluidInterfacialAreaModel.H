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
    Foam::fluidInterfacialAreaModel

Description
    Class to handle fluid-fluid and fluid-structure interfacial area.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    fluidInterfacialAreaModel.C
    newDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef fluidInterfacialAreaModel_H
#define fluidInterfacialAreaModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "stringOps.H"
#include "fluid.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidInterfacialAreaModel Declaration
\*---------------------------------------------------------------------------*/

class fluidInterfacialAreaModel
:
    public IOdictionary
{
protected:

    //- Reference to dispersed fluid
    const fluid& dispersed_;

    //- Reference to continuous fluid
    const fluid& continuous_;

    //- Adjusted pair residual
    const dimensionedScalar residualAlpha_;

public:

    //- Runtime type information
    TypeName("fluidInterfacialAreaModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidInterfacialAreaModel,
            fluidInterfacialAreaModels,
            (
                const objectRegistry& objReg,
                const dictionary& dict,
                const fluid& dispersed,
                const fluid& continuous
            ),
            (objReg, dict, dispersed, continuous)
        );

    //- Static data members
    static const dimensionSet dimIA;

    //- Constructors

        fluidInterfacialAreaModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous
        );

        //- Destructor
        virtual ~fluidInterfacialAreaModel();

    //- Selectors

        static autoPtr<fluidInterfacialAreaModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous
        );

    //- Member Functions

        //- Return interfacial area
        virtual tmp<volScalarField> iA() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
