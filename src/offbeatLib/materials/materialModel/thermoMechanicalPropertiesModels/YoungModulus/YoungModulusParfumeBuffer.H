/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::YoungModulusParfumeBuffer

Description
    Class modelling Young Modulus of Bufffer from Parfume.

    \htmlinclude YoungModulus_YoungModulusParfumeBuffer.html

SourceFiles
    YoungModulusParfumeBuffer.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef YoungModulusParfumeBuffer_H
#define YoungModulusParfumeBuffer_H

#include "YoungModulusModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class YoungModulusParfumeBuffer Declaration
\*---------------------------------------------------------------------------*/

class YoungModulusParfumeBuffer
:
    public YoungModulusModel
{
    // Private data

        //- Name of fast fluence field
        const word fastFluenceName_;

        //- Fast fluence field
        const volScalarField* phi_;

        //- Name of density field
        const word densityName_;

        //- Density field
        const volScalarField* rho_;

        //- Parameters for the Parfume model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar par6;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        YoungModulusParfumeBuffer(const YoungModulusParfumeBuffer&);

        //- Disallow default bitwise assignment
        void operator=(const YoungModulusParfumeBuffer&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("BufferParfume");

    // Constructors

        //- Construct from dictionary
        YoungModulusParfumeBuffer
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict,
            const word defaultModel
        );
    //- Destructor
    virtual ~YoungModulusParfumeBuffer();


    // Member Functions

    //- Update conductivity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
