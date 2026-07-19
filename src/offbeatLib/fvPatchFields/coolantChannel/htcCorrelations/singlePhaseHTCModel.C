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

#include "singlePhaseHTCModel.H"
#include "IF97.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

singlePhaseHTCModel::Model
singlePhaseHTCModel::parse(const word& name)
{
    if (name == "analyticalLaminar") return Model::analyticalLaminar;
    if (name == "Gnielinski")        return Model::Gnielinski;
    if (name == "DittusBoelter")     return Model::DittusBoelter;
    if (name == "TRACETubeNC")       return Model::TRACETubeNC;
    if (name == "TRACERodBundleNC")  return Model::TRACERodBundleNC;
    if (name == "ChurchillChu")      return Model::ChurchillChu;

    FatalErrorInFunction
        << "Unknown singlePhase correlation \"" << name << "\"." << nl
        << "Valid options: analyticalLaminar, Gnielinski, DittusBoelter, "
        << "TRACETubeNC, TRACERodBundleNC, ChurchillChu"
        << abort(FatalError);

    return Model::Gnielinski;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

singlePhaseHTCModel::singlePhaseHTCModel(const word& name)
:
    model_(parse(name))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar singlePhaseHTCModel::compute
(
    const scalar T0,
    const scalar Tw,
    const scalar p,
    const scalar G,
    const scalar DH
) const
{
    const scalar g = 9.80665;

    const scalar Tfilm = 0.5*(T0 + Tw);

    const scalar mu   = IF97::visc_Tp(Tfilm, p);
    const scalar Pr   = IF97::prandtl_Tp(Tfilm, p);
    const scalar k    = IF97::tcond_Tp(Tfilm, p);
    const scalar rho  = IF97::rhomass_Tp(Tfilm, p);
    const scalar beta = IF97::volumetricThermalExpansion_T(Tfilm);

    const scalar Re = G * DH / max(mu, SMALL);

    switch (model_)
    {
        case Model::analyticalLaminar:
        {
            return max(k/DH * 4.36, SMALL);
        }

        case Model::Gnielinski:
        {
            const scalar ReEff = max(Re, 2300.0);
            const scalar f = pow(1.58*log(ReEff) - 3.28, -2.0);

            return max
            (
                k/DH
              * ((0.5*f)*(ReEff - 1000.0)*Pr)
              / (1.0 + 12.7*sqrt(0.5*f)*(pow(Pr, 2.0/3.0) - 1.0)),
                SMALL
            );
        }

        case Model::DittusBoelter:
        {
            return max(k/DH * 0.023 * pow(Re, 0.8) * pow(Pr, 0.4), SMALL);
        }

        case Model::TRACETubeNC:
        {
            const scalar deltaT = mag(Tw - T0);
            const scalar nu = mu/max(rho, SMALL);

            const scalar hTurb =
                0.1 * k/DH * pow(Pr, 1.0/3.0)
              * cbrt(g * beta * deltaT * pow(DH, 3.0) / max(nu*nu, SMALL));

            const scalar hLam =
                0.59 * k/DH * pow(Pr, 0.25)
              * pow(g * beta * deltaT * pow(DH, 3.0) / max(nu*nu, SMALL), 0.25);

            return max(max(hTurb, hLam), SMALL);
        }

        case Model::TRACERodBundleNC:
        {
            const scalar deltaT = mag(Tw - T0);
            const scalar nu = mu/max(rho, SMALL);

            return max
            (
                0.7 * k/DH * pow(Pr, 0.25)
              * pow(g * beta * deltaT * pow(DH, 3.0) / max(nu*nu, SMALL), 0.25),
                SMALL
            );
        }

        case Model::ChurchillChu:
        {
            // Churchill & Chu (1975), valid for all Ra.
            // DH is the vertical plate height L.
            const scalar L = DH;
            const scalar deltaT = mag(Tw - T0);
            const scalar nu = mu/max(rho, SMALL);

            const scalar Ra =
                g * beta * deltaT * pow3(L) * Pr
              / max(sqr(nu), SMALL);

            const scalar num = 0.387 * pow(Ra, 1.0/6.0);
            const scalar den =
                pow(1.0 + pow(0.492/max(Pr, SMALL), 9.0/16.0), 8.0/27.0);

            return max(sqr(0.825 + num/den) * k / max(L, SMALL), SMALL);
        }
    }

    FatalErrorInFunction
        << "Unhandled singlePhase model."
        << abort(FatalError);

    return SMALL;
}

} // End namespace Foam

// ************************************************************************* //
