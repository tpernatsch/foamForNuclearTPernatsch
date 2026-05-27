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
    Foam::diffCoefPiecewiseArrhenius

Description
    Piecewise Arrhenius diffusion coefficient: uses parameters (d1, q1) below
    Tswitch and (d2, q2) at or above Tswitch.

    \f[
        D = \begin{cases}
            d_1 \exp(-q_1/RT) & T < T_\text{switch} \\
            d_2 \exp(-q_2/RT) & T \geq T_\text{switch}
        \end{cases}
    \f]

    Typical use: Kr in the UO2 kernel (Tswitch = 1773.15 K).

    \verbatim
    Kernel
    {
        type        piecewiseArrhenius;
        d1          1.06e-5; q1  272.7e3;
        d2          7.43e-3; q2  326.8e3;
        Tswitch     1773.15;
    }
    \endverbatim

SourceFiles
    diffCoefPiecewiseArrhenius.C

\*---------------------------------------------------------------------------*/

#ifndef diffCoefPiecewiseArrhenius_H
#define diffCoefPiecewiseArrhenius_H

#include "diffCoefModel.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class diffCoefPiecewiseArrhenius Declaration
\*---------------------------------------------------------------------------*/

class diffCoefPiecewiseArrhenius
:
    public diffCoefModel
{
    scalar d1_;
    scalar q1_;
    scalar d2_;
    scalar q2_;
    scalar Tswitch_;

    diffCoefPiecewiseArrhenius(const diffCoefPiecewiseArrhenius&) = delete;
    void operator=(const diffCoefPiecewiseArrhenius&) = delete;

public:

    TypeName("piecewiseArrhenius");

    diffCoefPiecewiseArrhenius
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const word defaultModel
    );

    virtual ~diffCoefPiecewiseArrhenius() {}

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
