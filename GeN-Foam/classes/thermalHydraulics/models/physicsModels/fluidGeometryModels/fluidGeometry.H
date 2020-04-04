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
#include "structureVolumetricAreaPartitionModel.H"

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

    const fluid& dispersed_;

    const fluid& continuous_;

    const structureModel& structure_;

    volScalarField DhDispersed_;

    volScalarField iA_;

    volScalarField frac1_;

    volScalarField frac2_;

    autoPtr<fluidDiameterModel> fluidDiameter_;

    autoPtr<fluidInterfacialAreaModel> fluidInterfacialArea_;

    autoPtr<structureVolumetricAreaPartitionModel>
        structureInterfacialAreaPartition_;

public:

    //- Constructors

        fluidGeometry
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const fluid& fluid1,
            const fluid& fluid2,
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

            const volScalarField& frac1() const
            {
                return frac1_;
            }

            const volScalarField& frac2() const
            {
                return frac2_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
