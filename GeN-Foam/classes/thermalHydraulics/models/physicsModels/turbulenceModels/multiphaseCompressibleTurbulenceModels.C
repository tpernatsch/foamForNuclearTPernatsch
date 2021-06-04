/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "phaseCompressibleTurbulenceModel.H"
#include "addToRunTimeSelectionTable.H"
#include "makeTurbulenceModel.H"

#include "laminarModel.H"
#include "RASModel.H"
#include "LESModel.H"
#include "rhoThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*
The command below creates, among others, the following typedef

typedef ThermalDiffusivity<PhaseCompressibleTurbulenceModel<rhoThermo>>
        rhoThermoPhaseCompressibleTurbulenceModel;

*/
makeTurbulenceModelTypes
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    rhoThermo
);

makeBaseTurbulenceModel
(
    volScalarField,
    volScalarField,
    compressibleTurbulenceModel,
    PhaseCompressibleTurbulenceModel,
    ThermalDiffusivity,
    rhoThermo
);

#define makeLaminarModel(Type)                                                 \
    makeTemplatedLaminarModel                                                  \
    (rhoThermoPhaseCompressibleTurbulenceModel, laminar, Type)

/*
The command below creates, among others, the following typedef

typedef Foam::RAS##Models::Type
        <
            Foam::EddyDiffusivity
            <   
                Foam::rhoThermoPhaseCompressibleTurbulenceModel
            >
        >
        Type##RAS##rhoThermoPhaseCompressibleTurbulenceModel;
*/
#define makeRASModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (rhoThermoPhaseCompressibleTurbulenceModel, RAS, Type)


#define makeLESModel(Type)                                                     \
    makeTemplatedTurbulenceModel                                               \
    (rhoThermoPhaseCompressibleTurbulenceModel, LES, Type)

#include "Stokes.H"
makeLaminarModel(Stokes);

#include "kEpsilon.H"
makeRASModel(kEpsilon);

//- In short, this creates the porousKEpsilon model with porousKEpsilon
//  templated with BasicTurbulenceModel = 
/*  
    EddyDiffusivity
    <
        ThermalDiffusivity
        <   
            PhaseCompressibleTurbulenceModel
            <
                rhoThermo
            >
        >
    >
*/
#include "porousKEpsilon.H"
makeRASModel(porousKEpsilon);

#include "LaheyKEpsilon.H"
makeRASModel(LaheyKEpsilon);

#include "mixtureKEpsilon.H"
makeRASModel(mixtureKEpsilon);


// ************************************************************************* //
