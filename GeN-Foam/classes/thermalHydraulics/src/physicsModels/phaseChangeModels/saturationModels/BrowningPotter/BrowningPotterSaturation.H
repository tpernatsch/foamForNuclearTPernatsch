/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::saturationModels::BrowningPotter

Description
    Saturation model based on results by Browning and Potter, refer to 
    https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf. The saturation pressure form is
    lnPSat = A + B/TSat + C*ln(TSat), which is not invertible in TSat. TSat is 
    inverted by approximateding ln(TSat) by a second order polyinomial.

SourceFiles
    BrowningPotterSaturation.C

\*---------------------------------------------------------------------------*/

#ifndef BrowningPotterSaturation_H
#define BrowningPotterSaturation_H

#include "saturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace saturationModels
{

/*---------------------------------------------------------------------------*\
                         Class BrowningPotter Declaration
\*---------------------------------------------------------------------------*/

class BrowningPotter
:
    public saturationModel
{
private :

    const volScalarField& iT_;

    const volScalarField& p_;

public:

    TypeName("BrowningPotter");

    //- Constructors

        //- Construct from a dictionary
        BrowningPotter
        (
            const phaseChangeModel& pcm,
            const dictionary& dict, 
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~BrowningPotter(){}

    //- Member Functions

        virtual scalar valueTSat(const label& celli) const;

        virtual scalar valuePSat(const label& celli) const;

        virtual scalar valueLnPSat(const label& celli) const;

        virtual scalar valuePSatPrime(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace saturationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
