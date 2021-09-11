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
    Foam::TONBModel

Description
    Model for the prediction of the temperature of the onset of nucleate
    boiling (i.e. T_ONB)

Author
    Stefan Radman (sradman@pm.me / stefan.radman@epfl.ch)

SourceFiles
    TONBModel.C
    newTONBModel.C

\*---------------------------------------------------------------------------*/

#ifndef TONBModel_H
#define TONBModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structure.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

//-
class FSPair;

/*---------------------------------------------------------------------------*\
                           Class TONBModel Declaration
\*---------------------------------------------------------------------------*/

class TONBModel
:
    public IOdictionary
{
protected:

    //- Ref to pair
    const FSPair& pair_;

    //- Ref to Tsat (i.e. the interface T)
    const volScalarField& Tsat_;

    //- Encapsulation sucks, doesn't it?
    //- Ptr to other fluid in pair
    mutable const fluid* otherFluidPtr_;

    //- Ptr to latent heat
    mutable const volScalarField* LPtr_;

public:

    TypeName("TONBModel");

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            TONBModel,
            TONBModels,
            (
                const FSPair& pair,
                const dictionary& dict,
                const objectRegistry& objReg
            ),
            (pair, dict, objReg)
        );

    //- Constructors

        TONBModel
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

        //- Destructor
        virtual ~TONBModel(){}

    //- Selectors

        static autoPtr<TONBModel> New
        (
            const FSPair& pair,
            const dictionary& dict,
            const objectRegistry& objReg
        );

    //- Member Functions

        //- Set ptrs
        void setPtrs() const;

        //- Return value of the temperature of the onset of nucleate boiling. 
        //  Unfortunately, this requires the 2-phase forced convection htc as
        //  well, and I have no way of getting that from registry
        virtual scalar value
        (
            const label& celli, 
            const scalar& htc2pFCi
        ) const = 0;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
