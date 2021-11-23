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

Class
    Foam::WallisFilm

Description
    Correlation to be used in annular flow scenarios in which the phase of
    interest to which this model is applied flows in a film along the 
    structures (e.g. pins) of interest. It assumes that the dimensionless
    film thickness d+ (=d/Dh with Dh being the structure hydraulic diameter)
    is proportional to the Wallis interfacial drag coefficient so that:
    
    0.005*75*alphaN = 0.005*300*d/Dh

    with alphaN being the normalized phase fraction of the phase of interest 
    and the film thickness is thus computed as:

    d = 0.25*alphaN*Dh

    See TRACE theory manual (https://www.nrc.gov/docs/ML1200/ML120060218.pdf,
    pp. 135-136, eqns. 4-37, 4-38) for further info  

Author
    Stefan Radman (stefan.radman@alumni.epfl.ch; sradman@protonmail.com)

SourceFiles
    WallisFilmFluidDiameter.C

\*---------------------------------------------------------------------------*/

#ifndef WallisFilm_H
#define WallisFilm_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidDiameterModel.H"

namespace Foam
{
namespace fluidDiameterModels
{

/*---------------------------------------------------------------------------*\
            Class WallisFilm Declaration
\*---------------------------------------------------------------------------*/

class WallisFilm
:
    public fluidDiameterModel
{
protected:

    //- Minimum allowable value of the phase fraction to avoid null fluid 
    //  diameter
    scalar residualAlpha_;
    
public:

    //- Runtime type information
    TypeName("WallisFilm");

    // Constructors

        WallisFilm
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~WallisFilm(){}

    // Member Functions

        //- Return fraction of structure interfacial contact in
        //  contact with the dispersed fluid
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fluidDiameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
