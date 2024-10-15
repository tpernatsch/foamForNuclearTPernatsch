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
    Foam::thermalExpansionPARFUMEBuffer

Description
    Class modelling the thermal expansion of Buffer through the code PARFUME.

    \htmlinclude thermalExpansion_thermalExpansionPARFUMEBuffer.html

SourceFiles
    thermalExpansionPARFUMEBuffer.C

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

#ifndef thermalExpansionPARFUMEBuffer_H
#define thermalExpansionPARFUMEBuffer_H

#include "thermalExpansionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermalExpansionPARFUMEBuffer Declaration
\*---------------------------------------------------------------------------*/

class thermalExpansionPARFUMEBuffer
:
    public thermalExpansionModel
{
    // Private data

        //- Parameters for the PARFUME model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thermalExpansionPARFUMEBuffer(const thermalExpansionPARFUMEBuffer&);

        //- Disallow default bitwise assignment
        void operator=(const thermalExpansionPARFUMEBuffer&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("BufferPARFUME");

    // Constructors

        //- Construct from dictionary
        thermalExpansionPARFUMEBuffer
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word defaultModel
        );

    //- Destructor
    virtual ~thermalExpansionPARFUMEBuffer();


    // Member Functions

    //- Update thermalExpansion
    virtual void correct(symmTensorField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
