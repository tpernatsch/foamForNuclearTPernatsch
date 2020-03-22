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
    Foam::dragModels::ChurchillGunterShaw

Description
    Churchill correlation for pressure drop along local principal axis 
    direction, Gunter-Shaw correlation for the pressure drop in the transverse 
    directions

SourceFiles
    ChurchillGunterShaw.C

\*---------------------------------------------------------------------------*/

#ifndef ChurchillGunterShaw_H
#define ChurchillGunterShaw_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class ChurchillGunterShaw Declaration
\*---------------------------------------------------------------------------*/

class ChurchillGunterShaw
:
    public dragModel
{
private:

    //- These three vector components tell me if the principal axis is aligned
    //  along the local X, Y or Z axes (local axis depend on the 
    //  structureModel). The principal component direction index is always
    //  stored on the last entry (i.e. [2]), while the transverse direction
    //  indices are stored in the first two entries. The direction indices
    //  are 0 for localX, 1 for localY, 2 for localZ
    vector lcmpt_;

    //- Surface roughness of the structure
    scalar surfaceRoughness_;
    
public:

    TypeName("ChurchillGunterShaw");

    //- Constructors

        //- Construct for a fluid-structure pair
        ChurchillGunterShaw
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~ChurchillGunterShaw(){}


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