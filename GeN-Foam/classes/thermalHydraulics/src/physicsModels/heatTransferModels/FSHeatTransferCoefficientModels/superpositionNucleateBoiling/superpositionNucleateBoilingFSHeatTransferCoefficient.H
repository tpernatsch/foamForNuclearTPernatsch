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
    Heat transfer coefficient for both single-phase convective scenarios that
    evolve into nucleate boiling ones, but no sub-cooled boiling region in
    between.

Dictionary entries

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    superpositionNucleateBoiling.C

\*---------------------------------------------------------------------------*/

#ifndef superpositionNucleateBoilingFSHeatTransferCoefficient_H
#define superpositionNucleateBoilingFSHeatTransferCoefficient_H

#include "FSHeatTransferCoefficientModel.H"
#include "flowEnhancementFactorModel.H"
#include "suppressionFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace FSHeatTransferCoefficientModels
{

/*---------------------------------------------------------------------------*\
                           Class superposition Declaration
\*---------------------------------------------------------------------------*/

class superpositionNucleateBoiling
:
    public FSHeatTransferCoefficientModel
{
private:

    //- Forced convection heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcFCPtr_;

    //- Pool boiling heat transfer coefficient
    autoPtr<FSHeatTransferCoefficientModel> htcPBPtr_;

    //- 
    autoPtr<flowEnhancementFactorModel> FPtr_;

    //- 
    autoPtr<suppressionFactorModel> SPtr_;

public:

    TypeName("superpositionNucleateBoiling");

    //- Constructors

        //- Construct for a fluid-structure pair
        superpositionNucleateBoiling
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~superpositionNucleateBoiling(){}


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