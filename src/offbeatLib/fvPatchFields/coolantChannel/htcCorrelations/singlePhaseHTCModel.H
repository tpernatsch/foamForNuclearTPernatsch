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

Class
    Foam::singlePhaseHTCModel

Description
    Single-phase forced- and natural-convection heat transfer correlations.

    Supported models:
    - analyticalLaminar    : Nu = 4.36 (constant wall heat flux, laminar)
    - Gnielinski           : Gnielinski (1976) for turbulent pipe flow
    - DittusBoelter        : Dittus-Boelter for turbulent pipe flow
    - TRACETubeNC          : TRACE natural-convection tube correlation
    - TRACERodBundleNC     : TRACE natural-convection rod-bundle correlation
    - ChurchillChu         : Churchill-Chu (1975) natural convection on a
                             vertical plate. DH is interpreted as the plate
                             height (vertical characteristic length).

SourceFiles
    singlePhaseHTCModel.C

\*---------------------------------------------------------------------------*/

#ifndef singlePhaseHTCModel_H
#define singlePhaseHTCModel_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class singlePhaseHTCModel Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseHTCModel
{
public:

    enum class Model
    {
        analyticalLaminar,
        Gnielinski,
        DittusBoelter,
        TRACETubeNC,
        TRACERodBundleNC,
        ChurchillChu
    };


    // Constructors

        singlePhaseHTCModel(const word& name);


    // Member Functions

        //- Compute HTC [W/(m²·K)].
        //  IF97 properties are evaluated at the film temperature
        //  0.5*(T0 + Tw) internally.
        scalar compute
        (
            scalar T0,   // bulk coolant temperature [K]
            scalar Tw,   // wall temperature         [K]
            scalar p,    // pressure                 [Pa]
            scalar G,    // mass flux kg/(m²·s)
            scalar DH    // hydraulic diameter       [m]
        ) const;

        Model model() const { return model_; }


private:

    Model model_;

    static Model parse(const word& name);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
