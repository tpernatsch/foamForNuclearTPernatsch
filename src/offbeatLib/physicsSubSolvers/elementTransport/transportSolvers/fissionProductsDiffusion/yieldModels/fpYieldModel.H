/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::fpYieldModel

Description
    Abstract base class for fission product yield (source term) models.

    Each subclass implements one production-rate correlation for a single FP
    species.  The production rate [mol/m³/s] is written into \c s for cells
    listed in \c addr.

    \verbatim
    Cs { yield { type enrichmentDependent; yieldLEU 0.14; yieldHEU 0.16; } }
    Kr { yield { type constant;            yield 0.297; } }
    Ag { yield { type burnupDependent;     yieldCoeff 1.31625e-3; burnupExponent 0.55734; } }
    \endverbatim

SourceFiles
    fpYieldModel.C

\*---------------------------------------------------------------------------*/

#ifndef fpYieldModel_H
#define fpYieldModel_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fpYieldModel Declaration
\*---------------------------------------------------------------------------*/

class fpYieldModel
{
    fpYieldModel(const fpYieldModel&) = delete;
    void operator=(const fpYieldModel&) = delete;

public:

    //- Runtime type information
    TypeName("fpYieldModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        fpYieldModel,
        dictionary,
        (const dictionary& dict, const scalar enrichment),
        (dict, enrichment)
    );

    // Constructors

        //- Construct from dictionary
        fpYieldModel(const dictionary& dict, const scalar enrichment)
        {}

    // Selector

        static autoPtr<fpYieldModel> New
        (
            const dictionary& dict,
            const scalar enrichment
        );

    //- Destructor
    virtual ~fpYieldModel()
    {}

    // Member Functions

        //- Write dimensionless yield Y [atoms/fission] for every cell.
        //  Solver applies s = Y * (Q/eJPerFission) / Na.
        virtual void updateYield
        (
            scalarField& Y,
            const scalarField& Bu
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
