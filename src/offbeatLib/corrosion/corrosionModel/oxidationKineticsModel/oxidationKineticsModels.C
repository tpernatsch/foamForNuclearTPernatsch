/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2025 OpenFOAM Foundation
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

#include "oxidationKineticsModels.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

#define addOxidationKineticsModel(type,name)                                \
defineTemplateTypeNameAndDebugWithName(type, name, 0);                      \
addToRunTimeSelectionTable(oxidationKineticsModel, type, dictionary);

// Define a new lowHighOxidationKineticsModel based on the supplied low and
// high temperature model classes
#define addLoHiOxidationKineticsModel(loType, hiType)                       \
typedef lowHighOxidationKineticsModel<loType, hiType>                       \
        loType##_##hiType##_oxidationKineticsModel;                         \
                                                                            \
addOxidationKineticsModel                                                   \
(                                                                           \
    loType##_##hiType##_oxidationKineticsModel,                             \
    (                                                                       \
        std::string(loType::typeName_()) + "|"                              \
      + std::string(hiType::typeName_())                                    \
    ).c_str()                                                               \
)

addLoHiOxidationKineticsModel(EPRI_KWU_CE, BakerJust);
addLoHiOxidationKineticsModel(EPRI_KWU_CE, LeistikovPraterCourtright);
addLoHiOxidationKineticsModel(EPRI_KWU_CE, CathcartPawelPraterCourtright);

addLoHiOxidationKineticsModel(EPRI_SLI, BakerJust);
addLoHiOxidationKineticsModel(EPRI_SLI, LeistikovPraterCourtright);
addLoHiOxidationKineticsModel(EPRI_SLI, CathcartPawelPraterCourtright);

addLoHiOxidationKineticsModel(MATPRO_CORROS_PWR, BakerJust);
addLoHiOxidationKineticsModel(MATPRO_CORROS_PWR, LeistikovPraterCourtright);
addLoHiOxidationKineticsModel(MATPRO_CORROS_PWR, CathcartPawelPraterCourtright);

addLoHiOxidationKineticsModel(MATPRO_CORROS_BWR, BakerJust);
addLoHiOxidationKineticsModel(MATPRO_CORROS_BWR, CathcartPawelPraterCourtright);
addLoHiOxidationKineticsModel(MATPRO_CORROS_BWR, LeistikovPraterCourtright);

} // End namespace Foam

// ************************************************************************* //
