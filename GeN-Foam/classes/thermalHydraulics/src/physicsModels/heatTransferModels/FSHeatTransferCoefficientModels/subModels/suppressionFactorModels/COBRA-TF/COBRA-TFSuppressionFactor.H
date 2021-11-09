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

Description

Author
    Stefan Radman (sradman@protonmail.com / stefan.radman@epfl.ch)

SourceFiles
    COBRA_TF.C

\*---------------------------------------------------------------------------*/

#ifndef COBRA_TFFlowEnhancementFactor_H
#define COBRA_TFFlowEnhancementFactor_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "suppressionFactorModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace suppressionFactorModels
{

/*---------------------------------------------------------------------------*\
            Class COBRA_TF Declaration
\*---------------------------------------------------------------------------*/

class COBRA_TF
:
    public suppressionFactorModel
{
protected:

    //- Ref to structure wall temperature
    const volScalarField& Twall_;

    //- Ref to fluid temperature
    const volScalarField& Tf_;

    //- Ref to saturation temperature
    const volScalarField& Tsat_;

public:

    //- Runtime type information
    TypeName("COBRA-TF");

    // Constructors

        COBRA_TF
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~COBRA_TF(){}

    // Member Functions

        //- Return fraction of structure interfacial contact in
        //  contact with the dispersed fluid
        virtual scalar value(const label& celli) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace suppressionFactorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
