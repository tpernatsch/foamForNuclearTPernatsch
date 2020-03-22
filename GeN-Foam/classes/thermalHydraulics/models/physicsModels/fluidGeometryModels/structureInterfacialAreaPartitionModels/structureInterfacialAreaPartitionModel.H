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
    Foam::structureInterfacialAreaPartitionModel

Description
    Class to handle fluid-fluid and fluid-structure interfacial area.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    structureInterfacialAreaPartitionModel.C
    newDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef structureInterfacialAreaPartitionModel_H
#define structureInterfacialAreaPartitionModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fluid.H"
#include "structureModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class structureInterfacialAreaPartitionModel Declaration
\*---------------------------------------------------------------------------*/

class structureInterfacialAreaPartitionModel
:
    public IOdictionary
{
protected:

    const fvMesh& mesh_;

    const fluid& fluid1_;

    const fluid& fluid2_;

    const structureModel& structure_;

public:

    //- Runtime type information
    TypeName("structureInterfacialAreaPartitionModel");

    //- Declare runtime construction
        
        declareRunTimeSelectionTable
        (
            autoPtr,
            structureInterfacialAreaPartitionModel,
            structureInterfacialAreaPartitionModels,
            (
                const objectRegistry& objReg,
                const dictionary& dict,
                const fluid& fluid1,
                const fluid& fluid2,
                const structureModel& structure
            ),
            (objReg, dict, fluid1, fluid2, structure)
        );

    // Constructors

        structureInterfacialAreaPartitionModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid1,
            const fluid& fluid2,
            const structureModel& structure
        );

        //- Destructor
        virtual ~structureInterfacialAreaPartitionModel();

    // Selectors

        static autoPtr<structureInterfacialAreaPartitionModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid1,
            const fluid& fluid2,
            const structureModel& structure
        );

    // Member Functions

        //- Return fraction of structure interfacial area in
        //  contact with fluid1
        virtual tmp<volScalarField> frac1() const = 0;

        //- Return fraction of structure interfacial area in
        //  contact with fluid2
        virtual tmp<volScalarField> frac2() const = 0;

        // Access
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
