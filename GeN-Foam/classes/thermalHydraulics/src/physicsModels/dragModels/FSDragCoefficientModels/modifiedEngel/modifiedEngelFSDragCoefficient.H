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
    Foam::FSDragCoefficientModels::modifiedEngel

Description
    modifiedEngel correlation for pressure drop

    @article{BUBELIS20083299,
    title = "Review and proposal for best fit of wire-wrapped fuel bundle 
        friction factor and pressure drop predictions using various existing 
        correlations",
    journal = "Nuclear Engineering and Design",
    volume = "238",
    number = "12",
    pages = "3299 - 3320",
    year = "2008",
    issn = "0029-5493",
    doi = "https://doi.org/10.1016/j.nucengdes.2008.06.024",
    url = "http://www.sciencedirect.com/science/article/pii/S0029549308003762",
    author = "E. Bubelis and M. Schikorr",
    abstract = "The aim of this paper is to provide an overview of the existing
        wire-wrapped fuel bundle friction factor/pressure drop correlations and
        to qualitatively evaluate which of the existing friction factor 
        correlations are the best in retracing the results of a large set of 
        the experimental data available on wire-wrapped fuel assemblies tested 
        under different coolant conditions."
    }

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    modifiedEngel.C

\*---------------------------------------------------------------------------*/

#ifndef modifiedEngelFSDragCoefficient_H
#define modifiedEngelFSDragCoefficient_H

#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FSPair;

namespace FSDragCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class modifiedEngel Declaration
\*---------------------------------------------------------------------------*/

class modifiedEngel
:
    public FSDragCoefficientModel
{
public:

    TypeName("modifiedEngel");

    //- Constructors

        //- Construct for a fluid-structure pair
        modifiedEngel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~modifiedEngel(){}

    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSDragCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //