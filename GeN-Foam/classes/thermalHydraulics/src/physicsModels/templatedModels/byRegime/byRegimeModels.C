/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "addToRunTimeSelectionTable.H"

#include "FFPair.H"
#include "FSPair.H"
#include "byRegimeModel.H"
#include "dispersionModel.H"
#include "fluidDiameterModel.H"
#include "interfacialAreaModel.H"
#include "contactPartitionModel.H"
#include "FFDragCoefficientModel.H"
#include "FSDragCoefficientModel.H"
#include "virtualMassCoefficientModel.H"
#include "twoPhaseDragMultiplierModel.H"
#include "FFHeatTransferCoefficientModel.H"
#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef byRegimeModel<scalar, dispersionModel> 
        byRegimeDispersion;

    addNamedToRunTimeSelectionTable
    (
        dispersionModel,
        byRegimeDispersion,
        dispersionModels,
        byRegime
    );

    typedef byRegimeModel<scalar, fluidDiameterModel> 
        byRegimeFluidDiameter;

    addNamedToRunTimeSelectionTable
    (
        fluidDiameterModel,
        byRegimeFluidDiameter,
        fluidDiameterModels,
        byRegime
    );

    typedef byRegimeModel<scalar, interfacialAreaModel> 
        byRegimeInterfacialArea;

    addNamedToRunTimeSelectionTable
    (
        interfacialAreaModel,
        byRegimeInterfacialArea,
        interfacialAreaModels,
        byRegime
    );

    typedef byRegimeModel<scalar, contactPartitionModel> 
        byRegimeContactPartition;

    addNamedToRunTimeSelectionTable
    (
        contactPartitionModel,
        byRegimeContactPartition,
        contactPartitionModels,
        byRegime
    );

    typedef byRegimeModel<scalar, FFDragCoefficientModel> 
        byRegimeFFDragCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FFDragCoefficientModel,
        byRegimeFFDragCoefficient,
        FFDragCoefficientModels,
        byRegime
    );

    typedef byRegimeModel<scalar, FSDragCoefficientModel> 
        byRegimeFSDragCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FSDragCoefficientModel,
        byRegimeFSDragCoefficient,
        FSDragCoefficientModels,
        byRegime
    );

    typedef byRegimeModel<scalar, virtualMassCoefficientModel> 
        byRegimeVirtualMassCoefficient;

    addNamedToRunTimeSelectionTable
    (
        virtualMassCoefficientModel,
        byRegimeVirtualMassCoefficient,
        virtualMassCoefficientModels,
        byRegime
    );

    typedef byRegimeModel<tensor, twoPhaseDragMultiplierModel> 
        byRegimeTwoPhaseDragMultiplier;

    addNamedToRunTimeSelectionTable
    (
        twoPhaseDragMultiplierModel,
        byRegimeTwoPhaseDragMultiplier,
        twoPhaseDragMultiplierModels,
        byRegime
    );

    typedef byRegimeModel<scalar, FFHeatTransferCoefficientModel> 
        byRegimeFFHeatTransferCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FFHeatTransferCoefficientModel,
        byRegimeFFHeatTransferCoefficient,
        FFHeatTransferCoefficientModels,
        byRegime
    );

    typedef byRegimeModel<scalar, FSHeatTransferCoefficientModel> 
        byRegimeFSHeatTransferCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel,
        byRegimeFSHeatTransferCoefficient,
        FSHeatTransferCoefficientModels,
        byRegime
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
