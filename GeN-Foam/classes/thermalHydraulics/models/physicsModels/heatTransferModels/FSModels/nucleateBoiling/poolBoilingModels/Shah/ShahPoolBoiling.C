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

#include "ShahPoolBoiling.H"
#include "addToRunTimeSelectionTable.H"
#include "fluid.H"
#include "structureModel.H"
#include "FSPair.H"
#include "heatTransferModel.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferModels
{
namespace nucleateBoilingSubModels
{
namespace poolBoilingModels
{
    defineTypeNameAndDebug(Shah, 0);
    addToRunTimeSelectionTable
    (
        poolBoilingModel, 
        Shah, 
        poolBoilingModels
    );
}
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransferModels::nucleateBoilingSubModels::poolBoilingModels::
Shah::Shah
(
    const heatTransferModel& htm,
    const objectRegistry& objReg,
    const FSPair& FSPair
)
:
    poolBoilingModel
    (
        htm,
        objReg,
        FSPair
    ),
    useExplicitHeatFlux_
    (
        this->lookupOrDefault<bool>("useExplicitHeatFlux", false)
    ),
    Twall_
    (
        (!useExplicitHeatFlux_) ?
        &FSPair.structure().Twall() :
        nullptr
    ),
    Tf_
    (
        (!useExplicitHeatFlux_) ?
        &FSPair.fluidRef().thermo().T() :
        nullptr
    ),
    q_
    (
        (useExplicitHeatFlux_) ? 
        &mesh_.lookupObject<volScalarField>("heatFlux.structure") :
        nullptr
    ),
    p_(mesh_.lookupObject<volScalarField>("p")),
    pCrit_(3.5e7), //- Specific to Sodium
    C0_(13.7),
    C1_(6.9),
    m0_(0.22),
    m1_(0.12),
    n_(0.7),
    exp_(1.0/(1.0-n_)),
    pR0_(5e-4),
    deltaPR_(1e-3)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

/*
IMPLEMENTATION NOTES

The original pool boiling heat transfer coefficient by Shah is the following 
form:

htcPB = C * q^n * pR^m

with C, n, m being constants, pR being the reduced pressure (i.e. the 
ration of the fluid pressure to the fluid pressure at its critical point)
and q being the heat flux between wall and fluid. Needless to say, this form
is not particularly suitable for a numerical implementation in its current
form as htcPB is required to compute q in the first place. In general:

q = htc*(Twall-Tf)

with htc being the total heat transfer coefficient between wall and fluid,
Twall being the wall temperature and Tf the fluid temperature. By assuming
that htc ~ htcPB (i.e. that most of the heat flux is due to the pool
boiling mechanism), one can tha equate:

htcPB = C * (htcPB*(Twall-Tf))^n * pR^m

which can be re-arranged to isolate htcPB on the left hand side as:

htcPB = (C * (Twall-Tf)^n * pR^m)^(1.0/(1.0-n))

Which is how the correlation is implemented here.
The assumption of htc ~ htcPB can be rather crude but it is widely employed
to yiled a feasable numerical implementation in other computer codes too.
The TRACE code is an example, even though the pool boiling correlation used
in TRACE is the Gorenflo correlation, not the Shah one. Nonetheless, it is
still a correlation in which htcPB is a function of the wall heat flux

A footnote on the parameters C, m, n. Shah provided these values:

C = 13.7, m = 0.22 if pR < 1e-3
C = 6.9, m = 0.12 if pR > 1e-3
n = 0.7 always

Clearly the discotinuity at pR = 1e-3 should be avoided for numerical
stabilitiy. For this reason, the values of C and m are linearly interpolated
in the range (pR0_, pR0_+deltaPR_). Currently, pR0_ = 5e-4 and deltaPR_ = 1e-3.
The choice of these values is arbitrary and should/could change in the future.

*/

void Foam::heatTransferModels::nucleateBoilingSubModels::poolBoilingModels::
Shah::correct()
{    
    scalar C;
    scalar m;

    if (!useExplicitHeatFlux_)
    {
        forAll(cellList_, i)
        {
            const label& celli(cellList_[i]);

            scalar deltaT((*Twall_)[celli]-(*Tf_)[celli]);
            
            if (deltaT > 0.0)
            {
                //- Reduced pressure
                scalar pR(p_[celli]/pCrit_);

                //- Interpolate C, m to avoid discontinuity at pR = 1e-3
                scalar f
                (
                    min
                    (
                        max
                        (
                            (pR-pR0_)/deltaPR_,
                            0
                        ),
                        1.0
                    )
                );

                C = f*C1_+(1.0-f)*C0_;
                m = f*m1_+(1.0-f)*m0_;

                htc_[i] = 
                (
                    pow
                    (
                        C*pow(pR, m)*pow(deltaT, n_),
                        exp_
                    )
                ); 
            }
            else
                htc_[i] = 0;
        }
    }
    else
    {
        forAll(cellList_, i)
        {
            const label& celli(cellList_[i]);

            const scalar& q((*q_)[celli]);
            
            if (q > 0.0)
            {
                //- Reduced pressure
                scalar pR(p_[celli]/pCrit_);

                //- Interpolate C, m to avoid discontinuity at pR = 1e-3
                scalar f
                (
                    min
                    (
                        max
                        (
                            (pR-pR0_)/deltaPR_,
                            0
                        ),
                        1.0
                    )
                );

                C = f*C1_+(1.0-f)*C0_;
                m = f*m1_+(1.0-f)*m0_;

                htc_[i] = 
                (
                    C*pow(q, n_)*pow(pR, m)
                ); 
            }
            else
                htc_[i] = 0;
        }
    }
}


// ************************************************************************* //
