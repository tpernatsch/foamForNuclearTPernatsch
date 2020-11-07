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
    Run-time selectable class to handle the suppression factor of nucleate
    boiling models.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    suppressionFactorModel.C

\*---------------------------------------------------------------------------*/

#ifndef suppressionFactorModel_H
#define suppressionFactorModel_H

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
                           Class suppressionFactorModel Declaration
\*---------------------------------------------------------------------------*/

class suppressionFactorModel
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

    //- List of cells in which this suppressionFactorModel exists
    const labelList& cellList_;

    //- Marker field for this power model, = 1 in all the cells that are
    //  listed in cellsList, and 0 elsewhere
    const scalarField& cellField_;

    //- The actual suppression factor
    scalarField suppressionFactor_;

public:

    //- Runtime type information
    TypeName("suppressionFactorModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        suppressionFactorModel,
        suppressionFactorModels,
        (
            const heatTransferModel& htm,
            const objectRegistry& objReg,
            const FSPair& FSPair
        ),
        (htm, objReg, FSPair)
    );

    // Constructors

        suppressionFactorModel
        (
            const heatTransferModel& htm,
            const objectRegistry& objReg,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~suppressionFactorModel(){}

    // Selectors

        static autoPtr<suppressionFactorModel> New
        (
            const heatTransferModel& htm,
            const objectRegistry& objReg,
            const FSPair& FSPair
        );

    // Members

        //- Update the suppressionFactor_
        virtual void correct() = 0;

        // Access

            const scalarField& suppressionFactor() const
            {
                return suppressionFactor_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleateBoiling
} // End namespace heatTransferModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
