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
    Foam::diffCoefSingleArrhenius

Description
    Single-term Arrhenius diffusion coefficient:

    \f[ D = d_1 \exp\!\left(-\frac{q_1}{RT}\right) \f]

    \verbatim
    Buffer { type singleArrhenius; d1 1e-8; q1 100e3; }
    \endverbatim

SourceFiles
    diffCoefSingleArrhenius.C

\*---------------------------------------------------------------------------*/

#ifndef diffCoefSingleArrhenius_H
#define diffCoefSingleArrhenius_H

#include "diffCoefModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class diffCoefSingleArrhenius Declaration
\*---------------------------------------------------------------------------*/

class diffCoefSingleArrhenius
:
    public diffCoefModel
{
    scalar d1_;
    scalar q1_;

    diffCoefSingleArrhenius(const diffCoefSingleArrhenius&) = delete;
    void operator=(const diffCoefSingleArrhenius&) = delete;

public:

    TypeName("singleArrhenius");

    diffCoefSingleArrhenius
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word defaultModel
    );

    virtual ~diffCoefSingleArrhenius() {}

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
