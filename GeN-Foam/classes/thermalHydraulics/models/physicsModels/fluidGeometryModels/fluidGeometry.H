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
    Foam::fluidGeometry

Description
    Class to handle geometric-related properties of the phase system.
    These properties consists in:
        -   Dh of the dispersed phase (Dh of continuous is assumed to be equal
            to Dh of structure/subStructures);
        -   interfacialArea between the two fluids
        -   interfacialArea repartition between each of the two fluids and the
            structure

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    fluidGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef fluidGeometry_H
#define fluidGeometry_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "fluid.H"
#include "structureModel.H"
#include "fluidDiameterModel.H"
#include "fluidInterfacialAreaModel.H"
#include "fluidPartitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidGeometry Declaration
\*---------------------------------------------------------------------------*/

class fluidGeometry
:
    public IOdictionary
{
protected:

    //- Ref to dispersed fluid
    const fluid& dispersed_;

    //- Ref to continuous fluid
    const fluid& continuous_;

    //- Ref to structure
    const structureModel& structure_;

    //- Characteristic dimension of the dispersed phase, computed by 
    //  fluidDiameter_
    volScalarField DhDispersed_;

    //- Fluid-fluid interfacial area, computed by fluidInterfacialArea_
    volScalarField iA_;

    //- Coefficients that weighs how much of the total two-phase drag is
    //  applied to the dispersed fluid, computed by dragPartition_
    volScalarField fKdD_;

    //- Coefficients that weighs how much of the total two-phase drag is
    //  applied to the continuous fluid, computed by dragPartition_
    volScalarField fKdC_;

    //- Coefficients that weighs how much of the total two-phase drag is
    //  applied to the dispersed fluid, computed by heatTransferPartition_
    volScalarField fHtcD_;

    //- Coefficients that weighs how much of the total two-phase drag is
    //  applied to the continuous fluid, computed by heatTransferPartition_
    volScalarField fHtcC_;

    //-
    autoPtr<fluidDiameterModel> fluidDiameter_;

    //- 
    autoPtr<fluidInterfacialAreaModel> fluidInterfacialArea_;

    //- This is set only if the same partitionModel is used for both drag
    //  and heat transfer
    autoPtr<fluidPartitionModel> fluidPartition_;

    //- Else these two are set
    autoPtr<fluidPartitionModel> dragPartition_;

    //- 
    autoPtr<fluidPartitionModel> heatTransferPartition_;

public:

    //- Constructors

        fluidGeometry
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& dispersed,
            const fluid& continuous,
            const structureModel& structure
        );

        //- Destructor
        virtual ~fluidGeometry(){}

    //- Static Member Functions

        static autoPtr<fluidGeometry> makeModel
        (
            const objectRegistry& srcObjReg,
            const objectRegistry& dstObjReg,
            const dictionary& dict
        );

    //- Member Functions

        //- Update interfacial area, dispersed phase diameter,
        //  structure interfacial area partition coeffs
        void correct();

        // Access

            //- Return fKd for the provided fluid name
            const volScalarField& fKd(const word& name) const;

            //- Return fHtc for the provided fluid name
            const volScalarField& fHtc(const word& name) const;

            const fluid& dispersed() const
            {
                return dispersed_;
            }

            const word& nameDispersed() const
            {
                return dispersed_.name();
            }

            const fluid& continuous() const
            {
                return continuous_;
            }

            const word& nameContinuous() const
            {
                return continuous_.name();
            }

            const volScalarField& DhStructure() const
            {
                return structure_.Dh();
            }

            const volScalarField& DhDispersed() const
            {
                return DhDispersed_;
            }

            const volScalarField& iA() const
            {
                return iA_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
