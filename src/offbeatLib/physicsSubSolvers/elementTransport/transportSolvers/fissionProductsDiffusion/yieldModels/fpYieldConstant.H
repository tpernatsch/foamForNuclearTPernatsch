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
    Foam::fpYieldConstant

Description
    Constant fission-product yield.

    \verbatim
    yield { type constant; yield 0.297; }
    \endverbatim

SourceFiles
    fpYieldConstant.C

\*---------------------------------------------------------------------------*/

#ifndef fpYieldConstant_H
#define fpYieldConstant_H

#include "fpYieldModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fpYieldConstant Declaration
\*---------------------------------------------------------------------------*/

class fpYieldConstant
:
    public fpYieldModel
{
    scalar yield_;   // atoms per fission

public:

    TypeName("constant");

    fpYieldConstant(const dictionary& dict, const scalar enrichment);
    ~fpYieldConstant() {}

    virtual void updateYield(scalarField& Y, const scalarField& Bu) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
