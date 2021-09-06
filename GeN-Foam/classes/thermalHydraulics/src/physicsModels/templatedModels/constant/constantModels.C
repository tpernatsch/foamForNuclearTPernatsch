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
#include "constantModel.H"
#include "latentHeatModel.H"
#include "fluidDiameterModel.H"
#include "interfacialAreaModel.H"
#include "contactPartitionModel.H"
#include "FFDragCoefficientModel.H"
#include "FSDragCoefficientModel.H"
#include "virtualMassCoefficientModel.H"
#include "FFHeatTransferCoefficientModel.H"
#include "FSHeatTransferCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef constantModel<scalar, latentHeatModel> 
        constantLatentHeat;

    addNamedToRunTimeSelectionTable
    (
        latentHeatModel,
        constantLatentHeat,
        latentHeatModels,
        constant
    );

    typedef constantModel<scalar, fluidDiameterModel> 
        constantFluidDiameter;

    addNamedToRunTimeSelectionTable
    (
        fluidDiameterModel,
        constantFluidDiameter,
        fluidDiameterModels,
        constant
    );

    typedef constantModel<scalar, interfacialAreaModel> 
        constantIntefacialArea;

    addNamedToRunTimeSelectionTable
    (
        interfacialAreaModel,
        constantIntefacialArea,
        interfacialAreaModels,
        constant
    );

    typedef constantModel<scalar, contactPartitionModel> 
        constantContactPartition;

    addNamedToRunTimeSelectionTable
    (
        contactPartitionModel,
        constantContactPartition,
        contactPartitionModels,
        constant
    );

    typedef constantModel<scalar, FFDragCoefficientModel> 
        constantFFDragCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FFDragCoefficientModel,
        constantFFDragCoefficient,
        FFDragCoefficientModels,
        constant
    );

    typedef constantModel<scalar, FSDragCoefficientModel> 
        constantFSDragCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FSDragCoefficientModel,
        constantFSDragCoefficient,
        FSDragCoefficientModels,
        constant
    );

    typedef constantModel<scalar, virtualMassCoefficientModel> 
        constantVirtualMassCoefficient;

    addNamedToRunTimeSelectionTable
    (
        virtualMassCoefficientModel,
        constantVirtualMassCoefficient,
        virtualMassCoefficientModels,
        constant
    );

    typedef constantModel<scalar, FFHeatTransferCoefficientModel> 
        constantFFHeatTransferCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FFHeatTransferCoefficientModel,
        constantFFHeatTransferCoefficient,
        FFHeatTransferCoefficientModels,
        constant
    );

    typedef constantModel<scalar, FSHeatTransferCoefficientModel> 
        constantFSHeatTransferCoefficient;

    addNamedToRunTimeSelectionTable
    (
        FSHeatTransferCoefficientModel,
        constantFSHeatTransferCoefficient,
        FSHeatTransferCoefficientModels,
        constant
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
