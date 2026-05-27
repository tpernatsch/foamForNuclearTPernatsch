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
    Foam::fpYieldEnrichmentDependent

Description
    Enrichment-dependent fission-product yield (e.g. Cs-137).

    Uses \c yieldLEU when enrichment < \c enrichmentThreshold,
    \c yieldHEU otherwise.

    \verbatim
    yield
    {
        type                 enrichmentDependent;
        yieldLEU             0.14;
        yieldHEU             0.16;
        enrichmentThreshold  0.175;   // optional, default 0.175
    }
    \endverbatim

SourceFiles
    fpYieldEnrichmentDependent.C

\*---------------------------------------------------------------------------*/

#ifndef fpYieldEnrichmentDependent_H
#define fpYieldEnrichmentDependent_H

#include "fpYieldModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class fpYieldEnrichmentDependent Declaration
\*---------------------------------------------------------------------------*/

class fpYieldEnrichmentDependent
:
    public fpYieldModel
{
    scalar yield_;   // resolved at construction from enrichment

public:

    TypeName("enrichmentDependent");

    fpYieldEnrichmentDependent(const dictionary& dict, const scalar enrichment);
    ~fpYieldEnrichmentDependent() {}

    virtual void updateYield(scalarField& Y, const scalarField& Bu) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
