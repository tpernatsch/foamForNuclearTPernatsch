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
    Foam::fpYieldBurnupDependent

Description
    Burnup-dependent fission-product yield following a power law (e.g. Ag-110m).

    \f[ y(b) = yieldCoeff \cdot b^{burnupExponent} \f]

    where \f$b\f$ is the local burnup in %FIMA, converted from the field
    \c Bu [MWd/tUO2] using:
    \f[ b = Bu \times 10^{-3} / (burnupConversion) \f]
    with default \c burnupConversion = 8.232 (= 0.881 × 9.3706, UO2).

    \verbatim
    yield
    {
        type              burnupDependent;
        yieldCoeff        1.31625e-3;
        burnupExponent    0.55734;
        burnupConversion  8.232;    // optional, default 8.232
    }
    \endverbatim

SourceFiles
    fpYieldBurnupDependent.C

\*---------------------------------------------------------------------------*/

#ifndef fpYieldBurnupDependent_H
#define fpYieldBurnupDependent_H

#include "fpYieldModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class fpYieldBurnupDependent Declaration
\*---------------------------------------------------------------------------*/

class fpYieldBurnupDependent
:
    public fpYieldModel
{
    scalar yieldCoeff_;
    scalar burnupExponent_;
    scalar burnupConversion_;   // MWd/tUO2 → %FIMA

public:

    TypeName("burnupDependent");

    fpYieldBurnupDependent(const dictionary& dict, const scalar enrichment);
    ~fpYieldBurnupDependent() {}

    virtual void updateYield(scalarField& Y, const scalarField& Bu) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
