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
    Foam::diffCoefDoubleArrhenius

Description
    Two-term Arrhenius diffusion coefficient:

    \f[ D = d_1 \exp\!\left(-\frac{q_1}{RT}\right)
          + d_2 \exp\!\left(-\frac{q_2}{RT}\right) \f]

    \verbatim
    Kernel { type doubleArrhenius; d1 5.6e-8; q1 209e3; d2 5.2e-4; q2 362e3; }
    \endverbatim

SourceFiles
    diffCoefDoubleArrhenius.C

\*---------------------------------------------------------------------------*/

#ifndef diffCoefDoubleArrhenius_H
#define diffCoefDoubleArrhenius_H

#include "diffCoefModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class diffCoefDoubleArrhenius Declaration
\*---------------------------------------------------------------------------*/

class diffCoefDoubleArrhenius
:
    public diffCoefModel
{
    scalar d1_;
    scalar q1_;
    scalar d2_;
    scalar q2_;

    diffCoefDoubleArrhenius(const diffCoefDoubleArrhenius&) = delete;
    void operator=(const diffCoefDoubleArrhenius&) = delete;

public:

    TypeName("doubleArrhenius");

    diffCoefDoubleArrhenius
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word defaultModel
    );

    virtual ~diffCoefDoubleArrhenius() {}

    virtual void updateCoef
    (
        volScalarField& sf,
        const volScalarField& T,
        const labelList& addr,
        const word FpName
    );
};

} // End namespace Foam

#endif

// ************************************************************************* //
