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

#include "leidenfrostModel.H"
#include "IF97.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

const scalar leidenfrostModel::g_ = 9.80665;


leidenfrostModel::Model
leidenfrostModel::parse(const word& name)
{
    if (name == "GroeneveldStewartCorrected")
        return Model::GroeneveldStewartCorrected;
    if (name == "GroeneveldStewart")
        return Model::GroeneveldStewart;
    if (name == "Ohnishi")   return Model::Ohnishi;
    if (name == "Gotovskij") return Model::Gotovskij;
    if (name == "Zuber")     return Model::Zuber;

    FatalErrorInFunction
        << "Unknown Leidenfrost correlation \"" << name << "\"." << nl
        << "Valid options: GroeneveldStewartCorrected, GroeneveldStewart, "
        << "Ohnishi, Gotovskij, Zuber"
        << abort(FatalError);

    return Model::GroeneveldStewartCorrected;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

leidenfrostModel::leidenfrostModel(const word& name)
:
    model_(parse(name))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar leidenfrostModel::compute
(
    const scalar T0,
    const scalar p,
    const scalar Q,
    const scalar hFB
) const
{
    const scalar Tsat = IF97::Tsat97(p);
    scalar TLeid = GREAT;

    switch (model_)
    {
        case Model::GroeneveldStewartCorrected:
        {
            TLeid =
                557.85
              + 44.1*p/1e6
              - 3.72*pow(p/1e6, 2.0)
              - Q*10000.0/(2.82 + 1.22*p/1e6);

            if (p > 9e6)
            {
                TLeid =
                    (
                        557.85
                      + 44.1*9.0
                      - 3.72*pow(9.0, 2.0)
                      - Q*10000.0/(2.82 + 1.22*9.0)
                      - IF97::Tsat97(9e6)
                    )
                  * (22064000.0 - p)/(22064000.0 - 9000000.0)
                  + Tsat;
            }

            return max(TLeid, Tsat + 0.01);
        }

        case Model::GroeneveldStewart:
        {
            TLeid = 557.85 + 44.1*p/1e6 - 3.72*pow(p/1e6, 2.0);

            if (p > 9e6)
            {
                TLeid =
                    (
                        557.85 + 44.1*9.0 - 3.72*pow(9.0, 2.0)
                      - IF97::Tsat97(9e6)
                    )
                  * (22064000.0 - p)/(22064000.0 - 9000000.0)
                  + Tsat;
            }

            return max(TLeid, Tsat + 0.01);
        }

        case Model::Ohnishi:
        {
            TLeid = Tsat + 350.0 + 5.1*(Tsat - T0);
            return max(TLeid, Tsat + 0.01);
        }

        case Model::Gotovskij:
        {
            TLeid = Tsat + 100.0 + 8.0*(Tsat - T0);
            return max(TLeid, Tsat + 0.01);
        }

        case Model::Zuber:
        {
            // Use fast approximation: T_Leid = T_sat + q_Leid / h_FB(T_wall).
            // When needsFilmBoilingTemperatureSolve() is true the BC instead
            // calls computeMinFlux() then solveFilmBoilingTemperature() for a
            // self-consistent result.
            TLeid = computeMinFlux(p)/max(hFB, SMALL) + Tsat;
            return max(TLeid, Tsat + 0.01);
        }
    }

    FatalErrorInFunction
        << "Unhandled Leidenfrost model."
        << abort(FatalError);

    return ROOTVGREAT;
}


scalar leidenfrostModel::computeMinFlux(const scalar p) const
{
    const scalar Tsat = IF97::Tsat97(p);

    return
        0.176
      * IF97::rhovap_p(p)
      * IF97::latentHeatVaporization_p(p)
      * pow
        (
            g_
          * IF97::sigma97(Tsat)
          * (IF97::rholiq_p(p) - IF97::rhovap_p(p))
          / pow(IF97::rholiq_p(p) + IF97::rhovap_p(p), 2.0),
            0.25
        );
}


bool leidenfrostModel::needsFilmBoilingTemperatureSolve() const
{
    return model_ == Model::Zuber;
}

} // End namespace Foam

// ************************************************************************* //
