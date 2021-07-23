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
    Foam::FSHeatTransferCoefficientModels::superposition

Description
    Heat transfer coefficient for that covers many regimes, from single phase
    convection, two-phase convection, sub-cooled and nucleate boiling, to
    film boiling (the latter, hopefully).

Dictionary entries

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    multiRegimeBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef multiRegimeBoilingFSHeatTransferCoefficient_H
#define multiRegimeBoilingFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"
#include "subCooledBoilingFractionModel.H"
#include "flowEnhancementFactorModel.H"
#include "suppressionFactorModel.H"
#include "TONBModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class superposition Declaration
\*---------------------------------------------------------------------------*/

class multiRegimeBoiling
:
    public FSHeatTransferCoefficientModel
{
private:

    //- Ref to wall temperature
    const volScalarField& Tw_;

    //- Ref to fluid temperature
    const volScalarField& Tf_;

    //- Ptr to FFPair (fUcK eNcApSuLaTi1o0N!11!1)
    mutable const FFPair* FFPairPtr_;

    //- Forced convection heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcFCPtr_;

    //- Pool boiling heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcPBPtr_;

    //- Flow enhcement factor model: increases the convective heat trasfer
    //  coefficient to account for bubble turbulence, increase in flow velocity
    //  due to the vapour generation
    autoPtr<flowEnhancementFactorModel> FPtr_;

    //- Suppression factor model: "turns-off" the pool boiling component
    //  as the boiling transition to film boiling
    autoPtr<suppressionFactorModel> SPtr_;

    //- Temperature of the (sub-cooled bulk) fluid for at which bubble
    //  detachment occurs
    autoPtr<subCooledBoilingFractionModel> SCBFPtr_;

    //- Temperature of the onset of nucleate boiling model
    autoPtr<TONBModel> TONBPtr_;

    //- Ptr to mass transfer term at the wall (of the phaseChange class) so
    //  that it can be set here to account e.g. for sub-cooled boiling
    mutable volScalarField* dmdtWPtr_;

    //-
    int dmdtSign_;

public:

    TypeName("multiRegimeBoiling");

    //- Constructors

        //- Construct for a fluid-structure pair
        multiRegimeBoiling
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~multiRegimeBoiling(){}


    //- Member Functions

        //- Drag coefficient
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace FSHeatTransferCoefficientModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //