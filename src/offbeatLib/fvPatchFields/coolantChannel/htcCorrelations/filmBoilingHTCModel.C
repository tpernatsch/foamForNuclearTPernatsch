/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "filmBoilingHTCModel.H"
#include "IF97.H"
#include "error.H"
#include "mathematicalConstants.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

const scalar filmBoilingHTCModel::g_ = 9.80665;


filmBoilingHTCModel::Model
filmBoilingHTCModel::parse(const word& name)
{
    if (name == "GroeneveldAnnular")  return Model::GroeneveldAnnular;
    if (name == "GroeneveldTubular")  return Model::GroeneveldTubular;
    if (name == "BishopSandbergTong") return Model::BishopSandbergTong;
    if (name == "Frederking")         return Model::Frederking;
    if (name == "Sakurai")            return Model::Sakurai;

    FatalErrorInFunction
        << "Unknown filmBoiling correlation \"" << name << "\"." << nl
        << "Valid options: GroeneveldAnnular, GroeneveldTubular, "
        << "BishopSandbergTong, Frederking, Sakurai"
        << abort(FatalError);

    return Model::Frederking;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

filmBoilingHTCModel::filmBoilingHTCModel(const word& name)
:
    model_(parse(name))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar filmBoilingHTCModel::compute
(
    const scalar T0,
    const scalar Tw,
    const scalar p,
    const scalar Q,
    const scalar G,
    const scalar DH
) const
{
    const scalar Tsat  = IF97::Tsat97(p);
    const scalar Tfilm = 0.5*(T0 + Tw);

    const scalar mu = IF97::visc_Tp(Tfilm, p);
    const scalar Pr = IF97::prandtl_Tp(Tfilm, p);

    const scalar Re = G * DH / max(mu, SMALL);

    scalar hFB = SMALL;

    switch (model_)
    {
        case Model::GroeneveldAnnular:
        case Model::GroeneveldTubular:
        {
            const scalar tempQ = max(Q, 0.0);

            scalar A, B, C, D;

            if (model_ == Model::GroeneveldAnnular)
            {
                A = 0.052; B = 0.688; C = 1.26; D = -1.06;
            }
            else
            {
                A = 0.00109; B = 0.989; C = 1.41; D = -1.15;
            }

            const scalar Y =
                max
                (
                    1.0
                  - 0.1 * pow
                    (
                        (1.0 - tempQ)*(IF97::rholiq_p(p)/IF97::rhovap_p(p) - 1.0),
                        0.4
                    ),
                    0.1
                );

            hFB =
                A * IF97::tcondvap_p(p) / DH
              * pow
                (
                    G * DH / IF97::viscvap_p(p)
                  * (tempQ + (1.0 - tempQ)*IF97::rhovap_p(p)/IF97::rholiq_p(p)),
                    B
                )
              * pow(Pr, C)
              * pow(Y, D);

            break;
        }

        case Model::BishopSandbergTong:
        {
            const scalar X =
                (Q <= 0.0 ? 0.0 :
                (Q >= 1.0 ? 1.0 :
                    1.0
                  /
                    (
                        1.0
                      + (1.0 - Q)/Q
                      * IF97::rhovap_p(p)/IF97::rholiq_p(p)
                    )));

            const scalar rhoEq =
                X*IF97::rhovap_p(p) + (1.0 - X)*IF97::rholiq_p(p);

            hFB =
                0.0193 * IF97::tcond_Tp(Tfilm, p) / DH
              * pow(Re, 0.8)
              * pow(Pr, 1.23)
              * pow(IF97::rhovap_p(p)/rhoEq, 0.68)
              * pow(IF97::rhovap_p(p)/IF97::rholiq_p(p), 0.068);

            break;
        }

        case Model::Frederking:
        {
            const scalar Tvap = 0.5*(Tw + Tsat);
            const scalar hfg =
                IF97::latentHeatVaporization_p(p)
              + 0.5*IF97::cpvap_p(p)*(Tw - Tsat);

            const scalar fsubc =
                1.0
              + 0.1
              * pow(IF97::rholiq_p(p)/IF97::rhovap_p(p), 0.75)
              * IF97::cpliq_p(p)
              * (Tsat - T0)
              / IF97::latentHeatVaporization_p(p);

            hFB =
                0.2 * fsubc
              * cbrt
                (
                    pow(IF97::tcond_Tp(Tvap, p), 2.0)
                  * hfg * g_ * IF97::rhovap_p(p)
                  * (IF97::rholiq_p(p) - IF97::rhovap_p(p))
                  / IF97::viscvap_p(p)
                  / max(Tw - Tsat, SMALL)
                );

            break;
        }

        case Model::Sakurai:
        {
            const scalar Tvap = 0.5*(Tw + Tsat);
            const scalar Tliq = 0.5*(T0 + Tsat);
            const scalar hcap = IF97::cpmass_Tp(Tvap, p)*(Tw - Tsat);
            const scalar Lz   = IF97::latentHeatVaporization_p(p) + 0.5*hcap;
            const scalar Rz   =
                sqrt
                (
                    (IF97::rhomass_Tp(Tvap, p)*IF97::visc_Tp(Tvap, p))
                  / (IF97::rhomass_Tp(Tliq, p)*IF97::visc_Tp(Tliq, p))
                );

            const scalar Prl = IF97::prandtl_Tp(Tliq, p);
            const scalar Prv = IF97::prandtl_Tp(Tvap, p);
            const scalar Sp  = hcap/(Lz*Prv);
            const scalar Sc  = IF97::cpmass_Tp(Tliq, p)*(Tsat - T0)/Lz;

            const scalar Az =
                pow(Sc, 3.0)/27.0
              + Rz*Rz*Sp*Prl*Sc/3.0
              + pow(Rz*Sp*Prl, 2.0)/4.0;

            const scalar Bz =
                -4.0/27.0*Sc*Sc
              + 2.0/3.0*Sp*Prl*Sc
              - 32.0/27.0*Sp*Prl*Rz*Rz
              + 0.25*pow(Sp*Prl, 2.0)
              + 2.0/27.0*pow(Sc, 3.0)/pow(Rz, 2.0);

            const scalar Cz = 0.5*Rz*Rz*Sp*Prl;

            const scalar Ez =
                cbrt(Az + Cz*sqrt(mag(Bz)))
              + cbrt(Az - Cz*sqrt(mag(Bz)))
              + Sc/3.0;

            const scalar lambdaTaylor =
                2.0*constant::mathematical::pi
              *
                sqrt
                (
                    mag
                    (
                        IF97::sigma97(Tsat)
                      /
                        (
                            g_
                          * max
                            (
                                SMALL,
                                IF97::rhomass_Tp(Tliq, p) - IF97::rhomass_Tp(Tvap, p)
                            )
                        )
                    )
                );

            const scalar Grv =
                g_
              * (IF97::rhomass_Tp(Tliq, p) - IF97::rhomass_Tp(Tvap, p))
              /
                (
                    IF97::rhomass_Tp(Tvap, p)
                  * pow(IF97::visc_Tp(Tvap, p)/IF97::rhomass_Tp(Tvap, p), 2.0)
                )
              * pow(lambdaTaylor, 3.0);

            const scalar Mz =
                (Grv*Prv*Lz*pow(Ez, 3.0))
              /
                (
                    pow(Rz*Prl*Sp, 2.0)
                  * hcap
                  * (1.0 + Ez/(Sp*Prl))
                );

            hFB =
                0.823 * IF97::tcond_Tp(Tvap, p) / DH
              * pow(mag(Mz), 0.25);

            break;
        }
    }

    return max(hFB, SMALL);
}

} // End namespace Foam

// ************************************************************************* //
