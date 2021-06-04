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
    Foam::dragModels::upgradedChengTodreas

Description
    upgradedChengTodreas correlation for BARE rod bundles (no wire wrapped).
    This is a very hasty implementation, just a stub for a quick test for a
    very demanding revisor...

    This is from Appendix A of Nuclear Engineering and Design vol. 335, pp.
    356-373 (2018)

SourceFiles
    upgradedChengTodreas.C

\*---------------------------------------------------------------------------*/

#ifndef upgradedChengTodreas_H
#define upgradedChengTodreas_H

#include "dragModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace dragModels
{

/*---------------------------------------------------------------------------*\
                           Class upgradedChengTodreas Declaration
\*---------------------------------------------------------------------------*/

class upgradedChengTodreas
:
    public dragModel
{
private:

    //- Pitch-diameter ratio
    scalar PbyD_;
    scalar aL_;
    scalar bL_;
    scalar cL_;
    scalar aT_;
    scalar bT_;
    scalar cT_;
    scalar ReL_;
    scalar ReT_;
    scalar logReTReL_;
    scalar CL_;
    scalar CT_;
    
public:

    TypeName("upgradedChengTodreas");

    //- Constructors

        //- Construct for a fluid-structure pair
        upgradedChengTodreas
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~upgradedChengTodreas(){}


    //- Member Functions

        //- Drag coefficient
        virtual void correctKd(volTensorField& Kd) const;

        //- Scalar drag factor fd in Kd = 0.5*fd*alpha*rho*magU/Dh
        scalar fd(const scalar& Re) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace dragModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //