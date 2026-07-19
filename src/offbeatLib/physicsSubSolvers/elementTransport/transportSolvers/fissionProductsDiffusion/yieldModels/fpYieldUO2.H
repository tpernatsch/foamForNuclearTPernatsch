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
    Foam::fpYieldUO2

Description
    Built-in fission-product yield model for a UO2 kernel.  Dispatches
    internally per species using standard IAEA HTGR benchmark defaults.

    \verbatim
    Cs  enrichmentDependent  yieldLEU 0.135  yieldHEU 0.16  threshold 0.175
    Ag  burnupDependent      yieldCoeff 1.31625e-3  burnupExponent 0.55734
    Kr  constant             yield 0.297
    Sr  constant             yield 0.10
    \endverbatim

    The species is identified by the \c fpName keyword, which the solver
    injects into the dict before calling \c New.

    \verbatim
    yield
    {
        Kernel { type UO2; }
    }
    \endverbatim

SourceFiles
    fpYieldUO2.C

\*---------------------------------------------------------------------------*/

#ifndef fpYieldUO2_H
#define fpYieldUO2_H

#include "fpYieldModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fpYieldUO2 Declaration
\*---------------------------------------------------------------------------*/

class fpYieldUO2
:
    public fpYieldModel
{
    word   fpName_;
    scalar yieldCs_;          // constant yield for Cs (enrichment-selected at ctor)
    scalar yieldAgCoeff_;     // burnup-dependent coefficient for Ag
    scalar yieldAgExp_;       // burnup exponent for Ag
    scalar yieldAgConv_;      // burnup conversion factor (MWd/tUO2 → %FIMA)
    scalar yieldKr_;          // constant yield for Kr
    scalar yieldSr_;          // constant yield for Sr

public:

    TypeName("UO2");

    fpYieldUO2(const dictionary& dict, const scalar enrichment);
    ~fpYieldUO2() {}

    virtual void updateYield(scalarField& Y, const scalarField& Bu) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
