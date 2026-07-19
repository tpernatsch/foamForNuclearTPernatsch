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

#include "onbModel.H"
#include "IF97.H"
#include "error.H"

namespace Foam
{

// Unit conversion constants used by Bergles-Rohsenow (USC units)
static const scalar PsiToPa     = 6894.757;
static const scalar BtuToJ      = 1055.06;
static const scalar FeetToMeters = 0.3048;


// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

onbModel::Model
onbModel::parse(const word& name)
{
    if (name == "Basu")            return Model::Basu;
    if (name == "BerglesRohsenow") return Model::BerglesRohsenow;

    FatalErrorInFunction
        << "Unknown ONB correlation \"" << name << "\"." << nl
        << "Valid options: Basu, BerglesRohsenow"
        << abort(FatalError);

    return Model::Basu;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

onbModel::onbModel(const word& name)
:
    model_(parse(name))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar onbModel::compute
(
    const scalar T0,
    const scalar p,
    const scalar kCoolant,
    const scalar hSP,
    const label  faceI
) const
{
    switch (model_)
    {
        case Model::Basu:
        {
            const scalar Tsat = IF97::Tsat97(p);

            const scalar phiAngle = 0.663;
            const scalar Fphi = 1.0 - exp(-pow(phiAngle, 3.0) - 0.5*phiAngle);

            const scalar dTONB =
                2.0*hSP*Tsat
              * IF97::sigma97(Tsat)
              /
                (
                    pow(Fphi, 2.0)
                  * IF97::rhovap_p(p)
                  * IF97::latentHeatVaporization_p(p)
                  * kCoolant
                );

            const scalar dTsub = Tsat - T0;

            return
                T0
              + 0.25
              * pow
                (
                    sqrt(dTONB) + sqrt(dTONB + 4.0*dTsub),
                    2.0
                );
        }

        case Model::BerglesRohsenow:
        {
            const scalar Tsat = IF97::Tsat97(p);

            scalar TONB = max(T0 + 1.0, Tsat + 1.0);

            const scalar a = 15.6 * pow(p/PsiToPa, 1.156);
            const scalar n = 2.3/pow(p/PsiToPa, 0.0234);
            const scalar unitConv = BtuToJ/3600.0/pow(FeetToMeters, 2.0);

            for (label iter = 0; iter < 50; ++iter)
            {
                const scalar dT    = max(TONB - Tsat, 1e-6);
                const scalar dTimp = dT*1.8;

                const scalar fx =
                    a * pow(dTimp, n) * unitConv
                  - hSP*(TONB - T0);

                const scalar dfx =
                    a * n * 1.8 * pow(dTimp, n - 1.0) * unitConv
                  - hSP;

                if (mag(fx) < 1e-6*max(hSP*(TONB - T0), 1.0))
                {
                    return TONB;
                }

                scalar delta = 0.0;

                if (mag(dfx) > SMALL)
                {
                    delta = -fx/dfx;
                }
                else
                {
                    delta = 0.5;
                }

                delta = max(min(delta, 20.0), -20.0);

                const scalar newTONB = max(TONB + delta, Tsat + 1e-6);

                if (mag(newTONB - TONB)/max(TONB, SMALL) < 1e-5)
                {
                    return newTONB;
                }

                TONB = newTONB;
            }

            WarningInFunction
                << "BerglesRohsenow ONB iteration did not converge on face "
                << faceI << ". Using TONB = " << TONB << endl;

            return TONB;
        }
    }

    FatalErrorInFunction
        << "Unhandled ONB model."
        << abort(FatalError);

    return GREAT;
}

} // End namespace Foam

// ************************************************************************* //
