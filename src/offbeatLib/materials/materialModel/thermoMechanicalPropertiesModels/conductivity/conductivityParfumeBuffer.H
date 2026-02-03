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
    Foam::conductivityParfumeBuffer

Description
    Buffer thermal conductivity model.
    
    \htmlinclude conductivity_conductivityParfumeBuffer.html

SourceFiles
    conductivityParfumeBuffer.C

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

#ifndef conductivityParfumeBuffer_H
#define conductivityParfumeBuffer_H

#include "conductivityModel.H"
#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityParfumeBuffer Declaration
\*---------------------------------------------------------------------------*/

class conductivityParfumeBuffer
:
    public conductivityModel
{
    // Private data

    //- Name used for the density field
    const word densityName_;

    //- Reference to density field
    const volScalarField* rho_;

    //- Initial density of Buffer
    scalar rho_init_;

    //- Theoretical density of Buffer
    scalar rho_theo_;

    //- Initial thermal conductivity of Buffer
    scalar k_init_;

    //- Thermal conductivity of Buffer at its theoretical density
    scalar k_theo_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    conductivityParfumeBuffer(const conductivityParfumeBuffer&);

    //- Disallow default bitwise assignment
    void operator=(const conductivityParfumeBuffer&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("BufferParfume");

    // Constructors

    //- Construct from dictionary
    conductivityParfumeBuffer
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const dictionary& baseDict,
        const word defaultModel
    );

    //- Destructor
    virtual ~conductivityParfumeBuffer();


    // Member Functions

    //- Update conductivity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
