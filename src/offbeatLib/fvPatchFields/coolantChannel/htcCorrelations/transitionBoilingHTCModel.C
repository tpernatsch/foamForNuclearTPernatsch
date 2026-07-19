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

#include "transitionBoilingHTCModel.H"
#include "IF97.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

transitionBoilingHTCModel::Model
transitionBoilingHTCModel::parse(const word& name)
{
    if (name == "McDonoughMilichKing") return Model::McDonoughMilichKing;
    if (name == "None")                return Model::None;

    FatalErrorInFunction
        << "Unknown transitionBoiling correlation \"" << name << "\"." << nl
        << "Valid options: McDonoughMilichKing, None"
        << abort(FatalError);

    return Model::McDonoughMilichKing;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

transitionBoilingHTCModel::transitionBoilingHTCModel(const word& name)
:
    model_(parse(name))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar transitionBoilingHTCModel::compute
(
    const scalar p,
    const scalar qWall,
    const scalar TCHF,
    const scalar qCHF,
    const scalar hFB
) const
{
    scalar hTB = hFB;

    switch (model_)
    {
        case Model::McDonoughMilichKing:
        {
            const scalar Tsat = IF97::Tsat97(p);

            hTB =
                qWall
              /
                (
                    (qCHF - qWall)/(4150.0*exp(3970000.0/p))
                  + TCHF - Tsat
                );

            break;
        }

        case Model::None:
        {
            hTB = hFB;
            break;
        }
    }

    return max(hTB, hFB);
}

} // End namespace Foam

// ************************************************************************* //
