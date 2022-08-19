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
    Foam::FSHeatTransferCoefficientModels::CachardLiquid

Description
    Heat transfer coefficient for Inverted Annular Flow - for the liquid phase. 
    This code cooperate with the CachardVapourFSHeatTransferCoefficient. 
    Warning, this model has been coded for rod bundle, especially for delta model,
    the thickness of the film of vapor. 

Author :
    Gauthier LAZARE - gauthier.lazare@gmail.com

SourceFiles
    CachardLiquidFSHeatTransferCoefficient.C

\*---------------------------------------------------------------------------*/

#ifndef CachardLiquidFSHeatTransferCoefficient_H
#define CachardLiquidFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class CachardLiquid Declaration
\*---------------------------------------------------------------------------*/

class CachardLiquid
:
    public FSHeatTransferCoefficientModel
{
private:

    //- Ref to structure wall temperature
    const volScalarField& Twall_;
 
    //- Ref to saturation temperature 
    const volScalarField& Tsat_;

    //- Ref to liquid temperature 
    const volScalarField& Tfluid_;

    //- Ref to normalized void fraction
    const volScalarField& alpha_;

    //- Ref to pressure
    const volScalarField& p_;

    //- Ref to normalized void fraction
    const volScalarField& rhoLiquid_;

    //-Ref to DR rod diameter
    const volScalarField&  rhoVapour_;

    //-Ref to kg kappa vapour
    const volScalarField&  kg_;

    //-Ref to mug mu vapour
    const volScalarField&  mug_;

    //-Ref to Dh vapour
    const volScalarField&  Dh_;
    
    //- Ref to the mass transfer rate - breaking encapsulation to reach MultiRegimeBoiling
    mutable volScalarField hGamma_;

    //- Reference parameter
    //- Stefan Boltzman Constant in W/m2/K
    scalar sigmaSB_;
    //- Gravity 
    scalar gravity_;

    //- Input parameters
    //- Wall emissivity 
    scalar epsilonWall_;
    //- Liquid emissivity
    scalar epsilonLiq_;

    

public:

    TypeName("CachardLiquid");

    //- Constructors

        //- Construct for a fluid-structure pair
        CachardLiquid
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~CachardLiquid(){}


    //- Member Functions

        //- Heat transfer coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //