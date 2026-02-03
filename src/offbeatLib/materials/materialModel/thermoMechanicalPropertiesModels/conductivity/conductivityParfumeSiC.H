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
    Foam::conductivityParfumeSiC

Description
    Conductivity of PyC material from Parfume manual

    \htmlinclude conductivity_conductivityParfumeSiC.html



SourceFiles
    conductivityParfumeSiC.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    December 2022

\*---------------------------------------------------------------------------*/

#ifndef conductivityParfumeSiC_H
#define conductivityParfumeSiC_H

#include "conductivityModel.H"
#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"

#include "userParameters.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class conductivityParfumeSiC Declaration
\*---------------------------------------------------------------------------*/

class conductivityParfumeSiC
:
    public conductivityModel
{
    // Private data

    //- Parameters for the Parfume model
    scalar par1;
    scalar par2;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    conductivityParfumeSiC(const conductivityParfumeSiC&);

    //- Disallow default bitwise assignment
    void operator=(const conductivityParfumeSiC&);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("SiCParfume");

    // Constructors

    //- Construct from dictionary
    conductivityParfumeSiC
    (
        const fvMesh& mesh,
        const dictionary& dict,
        const dictionary& baseDict,
        const word defaultModel
    );

    //- Destructor
    virtual ~conductivityParfumeSiC();


    // Member Functions

    //- Update conductivity
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
