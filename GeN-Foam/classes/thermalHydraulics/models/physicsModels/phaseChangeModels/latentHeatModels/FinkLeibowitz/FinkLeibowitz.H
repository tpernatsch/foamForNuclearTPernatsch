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
    Foam::latentHeatModels::FinkLeibowitz

Description
    This class describes the latent heat of vaporization of liquid sodium
    according to the data presented in the report by Fink & Leibowitz 
    (ANL RE 95/2), https://www.ne.anl.gov/eda/ANL-RE-95-2.pdf

SourceFiles
    FinkLeibowitz.C

\*---------------------------------------------------------------------------*/

#ifndef FinkLeibowitz_H
#define FinkLeibowitz_H

#include "latentHeatModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace latentHeatModels
{

/*---------------------------------------------------------------------------*\
            Class FinkLeibowitz Declaration
\*---------------------------------------------------------------------------*/

class FinkLeibowitz
:
    public latentHeatModel
{
protected:


public:

    TypeName("FinkLeibowitz");

    //- Constructors

        //- Construct from a dictionary
        FinkLeibowitz
        (
            const dictionary& dict, 
            const fluid& fluid1,
            const fluid& fluid2,
            const volScalarField& p,
            const volScalarField& dmdt,
            const volScalarField& iT
        );

        //- Destructor
        virtual ~FinkLeibowitz(){}

    //- Member Functions

        virtual void correct();

        void setL();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace latentHeatModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
