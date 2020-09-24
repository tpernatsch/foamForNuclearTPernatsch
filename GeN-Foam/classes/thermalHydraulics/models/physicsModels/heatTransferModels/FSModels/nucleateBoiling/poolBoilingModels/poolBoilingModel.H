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
    Foam::poolBoilingModel

Description
    Run-time selectable class to handle the heat transfer coefficient
    associated with pure pool boiling (i.e. boiling without forced convection)

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    poolBoilingModel.C

\*---------------------------------------------------------------------------*/

#ifndef poolBoilingModel_H
#define poolBoilingModel_H

#include "dictionary.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;
class structureModel;
class FSPair;
class heatTransferModel;

namespace heatTransferModels
{

namespace nucleateBoilingSubModels
{

/*---------------------------------------------------------------------------*\
                           Class poolBoilingModel Declaration
\*---------------------------------------------------------------------------*/

class poolBoilingModel
:
    public IOdictionary
{
protected:

    //- Reference to mesh
    const fvMesh& mesh_;

    //- Local object registry (i.e. the regime to which the parent 
    //  heatTransferModel belongs to)
    const objectRegistry& objReg_;

    //- The fluid-structure pair of the parent heatTransferModel
    const FSPair& FSPair_;

    //- List of cells in which this poolBoilingModel exists
    const labelList& cellList_;

    //- Marker field for this power model, = 1 in all the cells that are
    //  listed in cellsList, and 0 elsewhere
    const scalarField& cellField_;

    //- Pool boiling heat transfer coefficient
    scalarField htc_;

public:

    //- Runtime type information
    TypeName("poolBoilingModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        poolBoilingModel,
        poolBoilingModels,
        (
            const heatTransferModel& htm,
            const objectRegistry& objReg,
            const FSPair& FSPair
        ),
        (htm, objReg, FSPair)
    );

    // Constructors

        poolBoilingModel
        (
            const heatTransferModel& htm,
            const objectRegistry& objReg,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~poolBoilingModel(){}

    // Selectors

        static autoPtr<poolBoilingModel> New
        (
            const heatTransferModel& htm,
            const objectRegistry& objReg,
            const FSPair& FSPair
        );

    // Members

        //- Update the poolBoiling_
        virtual void correct() = 0;

        // Access

            const scalarField& htc() const
            {
                return htc_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleateBoiling
} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
