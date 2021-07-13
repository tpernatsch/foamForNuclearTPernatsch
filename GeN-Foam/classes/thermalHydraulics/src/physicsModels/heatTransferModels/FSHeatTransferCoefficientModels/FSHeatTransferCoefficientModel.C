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

\*---------------------------------------------------------------------------*/

#include "FSPair.H"
#include "FSHeatTransferCoefficientModel.H"
#include "myOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FSHeatTransferCoefficientModel, 0);
    defineRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel, 
        FSHeatTransferCoefficientModels
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FSHeatTransferCoefficientModel::FSHeatTransferCoefficientModel
(
    const FSPair& pair,
    const dictionary& dict,
    const objectRegistry& objReg
)
:
    IOdictionary
    (
        IOobject
        (
            "FSHeatTransferCoefficientModel."+typeName,
            pair.mesh().time().timeName(),
            objReg,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dict
    ),
    pair_(pair),
    cellsPtr_(nullptr)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FSHeatTransferCoefficientModel::correctField(volScalarField& htc)
{
    if (cellsPtr_ == nullptr)
    {
        cellsPtr_ = new DynamicList<label>;
        DynamicList<label>& cells(*cellsPtr_);
        wordList regions(myOps::split(this->dictName(), ':'));
        forAll(regions, i)
        {
            const labelList& regionCells
            (
                pair_.structureRef().cellLists()[regions[i]]
            );
            forAll(regionCells, j)
            {
                cells.append(regionCells[j]);
            }
        }
    }
    const DynamicList<label>& cells(*cellsPtr_);
    forAll(cells, i)
    {
        const label& celli(cells[i]);
        htc[celli] = this->value(celli);
    }
    htc.correctBoundaryConditions();
}

// ************************************************************************* //
