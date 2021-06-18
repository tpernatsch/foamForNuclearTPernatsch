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

Description
    Base class to represent a FSDragFactor with basing name access
    functionality.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::FSDragFactor

SourceFiles
    FSDragFactor.C

\*---------------------------------------------------------------------------*/

#ifndef FSDragFactor_H
#define FSDragFactor_H

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "HashPtrTable.H"
#include "FSDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FSPair;

/*---------------------------------------------------------------------------*\
                           Class FSDragFactor Declaration
\*---------------------------------------------------------------------------*/

class FSDragFactor
    : public IOdictionary
{
protected :

    //- Ref to FSPair
    const FSPair& pair_;

    //- Cells in which this model operates
    labelList cells_;

    //- True if the drag factor tensor is isotropic 
    bool isotropic_;

    //- Structure phase fraction
    const volScalarField& alphas_;

    //- Normalized fluid phase fraction
    const volScalarField& alphaN_;

    //- Fluid density
    const volScalarField& rho_;

    //- Hydraulic diameter/fluid characteristic dimension
    const volScalarField& Dh_;

    //- Fluid velocity
    const volVectorField& U_;
    const volScalarField& magU_;

    //-
    mutable scalarField halfAlphaRhoMagUByDh_;
    
    //- Drag coefficient for an isotropiuc drag factor
    autoPtr<FSDragCoefficientModel> fdPtr_;

    //- Table of drag coefficients for an anisotropic drag factor
    typedef HashTable
    <
        autoPtr<FSDragCoefficientModel>,
        word,
        word::hash
    > fdTable;
    fdTable fdPtrs_;

    //-
    bool validX_;
    bool validY_;
    bool validZ_;

public:

    TypeName("FSDragFactor");

    //- Constructors

        FSDragFactor
        (
            const FSPair& pair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~FSDragFactor() {}

    //- Member Functions

        //- Drag coefficient used in momentum equations:
        //    ddt(alpha_i*rho_i*U_i) + ... = ... -Kd&U_i
        //  for fluid-structure drag.
        void correctField(volTensorField& Kd) const;

        //- Access

            /*
            const labelList& cells() const
            {
                return cells_;
            }*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
