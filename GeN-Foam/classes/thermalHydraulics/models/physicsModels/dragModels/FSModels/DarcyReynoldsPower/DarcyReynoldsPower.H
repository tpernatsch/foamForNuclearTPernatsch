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
    Foam::dragModels::DarcyReynoldsPower

Description
    Drag coefficient in the form coeff*Re^exp. The coeff and exp are specified
    in vector form, one per local component (i.e. in the local structure frame
    of reference). An isotropic variant can be specified via the isotropic
    flag, if the isotropic Reynolds and hydrauluic diameters are to be used
    to compute a single drag coefficient for all directions

SourceFiles
    DarcyReynoldsPower.C

\*---------------------------------------------------------------------------*/

#ifndef DarcyReynoldsPower_H
#define DarcyReynoldsPower_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class DarcyReynoldsPower Declaration
\*---------------------------------------------------------------------------*/

class DarcyReynoldsPower
:
    public dragModel
{
private:

    //- Coefficients in f = coeffs_.*Re.^(exps_)
    vector coeff_;

    //- Exponents in f = coeffs_.*Re.^(exps_)
    vector exp_;

    //- True if coeff_, exp_ uniform in all directions
    bool isotropic_;

public:

    TypeName("DarcyReynoldsPower");

    //- Constructors

        //- Construct for a fluid-structure pair
        DarcyReynoldsPower
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~DarcyReynoldsPower(){}


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