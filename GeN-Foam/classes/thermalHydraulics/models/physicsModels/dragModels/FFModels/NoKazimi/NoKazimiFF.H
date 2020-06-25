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
    Foam::dragModels::NoKazimi

Description
    NoKazimi correlation for interfacial liquid-vapour sodium friction
    https://dspace.mit.edu/handle/1721.1/60581. It was devisied for
    liquid-vapour sodium flows in pin bundles, and requires the
    specification of the bundle pin pitch and diameter
    
Dictionary entries (mandatory)
    pinPitch        0.006;
    pinDiameter     0.0076;
    
Author
    By Stefan Radman (stefan.radman@epfl.ch)

Source files
    NoKazimiFF.C

\*---------------------------------------------------------------------------*/

#ifndef NoKazimiFF_H
#define NoKazimiFF_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class NoKazimiFF Declaration
\*---------------------------------------------------------------------------*/

class NoKazimiFF
:
    public dragModel
{
protected:

    //- Ref to vapour
    const fluid& vapour_;

    //- Ref to liquid
    const fluid& liquid_;

    //- Pre-computed coeffs
    scalar A_;

public:

    TypeName("NoKazimi");

    //- Constructors

        //- Construct from a dictionary and a phase pair
        NoKazimiFF
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair
        );

        //- Destructor
        virtual ~NoKazimiFF(){}

    //- Member Functions

        //- Drag coefficient
        virtual void correctKd(volTensorField& Kd) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
