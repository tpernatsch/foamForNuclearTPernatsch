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

#include "saturatedBoilingHTCModel.H"
#include "IF97.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

saturatedBoilingHTCModel::Model
saturatedBoilingHTCModel::parse(const word& name)
{
    if (name == "Chen")            return Model::Chen;
    if (name == "SchrockGrossman") return Model::SchrockGrossman;
    if (name == "Rohsenow")        return Model::Rohsenow;

    FatalErrorInFunction
        << "Unknown saturatedBoiling correlation \"" << name << "\"." << nl
        << "Valid options: Chen, SchrockGrossman, Rohsenow"
        << abort(FatalError);

    return Model::Rohsenow;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

saturatedBoilingHTCModel::saturatedBoilingHTCModel(const word& name)
:
    model_(parse(name))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar saturatedBoilingHTCModel::compute
(
    const scalar T0,
    const scalar Tw,
    const scalar p,
    const scalar Q,
    const scalar G,
    const scalar DH,
    const scalar qWall,
    const scalar hSP
) const
{
    const scalar g = 9.80665;
    const scalar Tsat = IF97::Tsat97(p);
    const scalar Tfilm = 0.5*(T0 + Tw);
    const scalar tempQ = max(Q, 0.0);

    scalar hSB = hSP;

    switch (model_)
    {
        case Model::Chen:
        {
            const scalar recXtt =
                pow(tempQ/max(1.0 - tempQ, SMALL), 0.9)
              * pow(IF97::rholiq_p(p)/IF97::rhovap_p(p), 0.5)
              * pow(IF97::viscvap_p(p)/IF97::viscliq_p(p), 0.1);

            const scalar hc =
                0.023
              * pow((G * DH * (1.0 - tempQ)) / IF97::viscliq_p(p), 0.8)
              * pow(IF97::prandtlliq_p(p), 0.4)
              * IF97::tcondliq_p(p) / DH;

            scalar F = 1.0;
            if (recXtt >= 0.1)
            {
                F = 2.35 * pow(0.213 + recXtt, 0.736);
            }

            const scalar S =
                1.0
              /
                (
                    1.0
                  + 2.53e-6
                  * pow
                    (
                        G * DH / IF97::visc_Tp(Tfilm, p)
                      * (1.0 - tempQ) * pow(F, 1.25),
                        1.17
                    )
                );

            const scalar hNB =
                0.00122
              * pow(max(mag(Tw - Tsat), SMALL), 0.24)
              * pow(max(mag(IF97::psat97(Tw) - IF97::psat97(Tsat)), SMALL), 0.75)
              *
                (
                    pow(IF97::tcondliq_p(p), 0.79)
                  * pow(IF97::cpliq_p(p), 0.45)
                  * pow(IF97::rholiq_p(p), 0.49)
                )
              /
                (
                    pow(IF97::sigma97(Tsat), 0.5)
                  * pow(IF97::viscliq_p(p), 0.29)
                  * pow(IF97::latentHeatVaporization_p(p), 0.24)
                  * pow(IF97::rhovap_p(p), 0.24)
                );

            hSB = F*hc + S*hNB;
            break;
        }

        case Model::SchrockGrossman:
        {
            const scalar recXtt =
                pow(tempQ/max(1.0 - tempQ, SMALL), 0.9)
              * pow(IF97::rholiq_p(p)/IF97::rhovap_p(p), 0.5)
              * pow(IF97::viscvap_p(p)/IF97::viscliq_p(p), 0.1);

            hSB =
                (
                    7400.0 * qWall / IF97::latentHeatVaporization_p(p) / G
                  + 1.11 * pow(recXtt, 0.66)
                )
              * hSP;
            break;
        }

        case Model::Rohsenow:
        {
            hSB =
                qWall
              /
                (
                    cbrt
                    (
                        qWall
                      /
                        (
                            IF97::viscliq_p(p)
                          * IF97::latentHeatVaporization_p(p)
                          * pow
                            (
                                g
                              * (IF97::rholiq_p(p) - IF97::rhovap_p(p))
                              / IF97::sigma97(Tsat),
                                0.5
                            )
                        )
                    )
                  *
                    (
                        0.013
                      * IF97::latentHeatVaporization_p(p)
                      * pow(IF97::prandtlliq_p(p), 1.0)
                    )
                  / IF97::cpliq_p(p)
                  + Tsat - T0
                );
            break;
        }
    }

    return max(hSB, hSP);
}

} // End namespace Foam

// ************************************************************************* //
