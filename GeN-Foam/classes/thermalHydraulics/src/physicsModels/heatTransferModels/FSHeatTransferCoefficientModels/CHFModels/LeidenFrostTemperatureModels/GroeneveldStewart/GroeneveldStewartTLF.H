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
    Foam::TLFModel::GroeneveldStewartTLF

Description
    Calculate the Leidenfrost Temperature at which the heat flux is minimal in the
    hysteresis after the CHF. The formula used is the same as TRACE with the correlation 
    of GroeneveldStewart. WARNING, it is not complete yet, it lack the term due to
    quality.

Author
    Gauthier Lazare (gauthier.lazare@gmail.com)

SourceFiles
    GroeneveldStewartTLF.C

\*---------------------------------------------------------------------------*/

#ifndef GroeneveldStewartTLF_H
#define GroeneveldStewartTLF_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TLFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace LeidenFrostTemperatureModels
{

/*---------------------------------------------------------------------------*\
            Class GroeneveldStewart Declaration
\*---------------------------------------------------------------------------*/

class GroeneveldStewart
:
    public TLFModel
{
protected:

    //- Critical Pressure of the liquid (Pa)
    scalar criticalPressure_;

public:

    //- Runtime type information
    TypeName("GroeneveldStewart");

    // Constructors

        GroeneveldStewart
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~GroeneveldStewart(){}

    // Member Functions

        //- Return TLF
        virtual scalar value
        (
            const label& celli
        ) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace LeidenFrostTemperatureModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
