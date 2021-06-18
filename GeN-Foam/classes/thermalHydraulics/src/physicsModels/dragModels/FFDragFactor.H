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
    Base class to represent a FFDragFactor with basing name access
    functionality.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::FFDragFactor

SourceFiles
    FFDragFactor.C

\*---------------------------------------------------------------------------*/

#ifndef FFDragFactor_H
#define FFDragFactor_H

#include "volFields.H"
#include "zeroGradientFvPatchFields.H"
#include "HashPtrTable.H"
#include "FFDragCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
                           Class FFDragFactor Declaration
\*---------------------------------------------------------------------------*/

class FFDragFactor
    : public IOdictionary
{
protected :

    //- Ref to FFPair
    const FFPair& pair_;

    //- Registry where FFDragCoefficients are stored to avoid confliciting
    //  FFDragCoefficient IOobject names if storing them all in the mesh
    //  registry (is it even a problem? Either way, it's cleaner)
    objectRegistry coeffReg_;

    //- Continutous phase density
    const volScalarField& rhoc_;

    //- Characteristic dimension of dispersed phase
    const volScalarField& Dhd_;

    //- Slip velocity
    const volScalarField& magUr_;
    
    //- Drag coefficient for an isotropiuc drag factor
    autoPtr<FFDragCoefficientModel> fdPtr_;

public:

    TypeName("FFDragFactor");

    //- Constructors

        FFDragFactor
        (
            const FFPair& pair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~FFDragFactor() {}

    //- Member Functions

        //- Drag coefficient used in momentum equations:
        //    ddt(alpha_i*rho_i*U_i) + ... = ... -Kd&U_i
        //  for fluid-structure drag.
        void correctField(volScalarField& Kd) const;

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
