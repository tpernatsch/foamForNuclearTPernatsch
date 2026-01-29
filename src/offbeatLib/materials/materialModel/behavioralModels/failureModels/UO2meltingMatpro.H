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
    Foam::UO2meltingMatpro

Description
    Model to check failure of UO2 based on a limit melting temperature from
    Matprov11.

Usage


SourceFiles
    UO2meltingMatpro.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef UO2meltingMatpro_H
#define UO2meltingMatpro_H

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class UO2meltingMatpro Declaration
\*---------------------------------------------------------------------------*/

class UO2meltingMatpro
:
    public failureModel
{
    // Private data

        //- Melting temperature UO2
        const scalar Tmelt_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        UO2meltingMatpro(const UO2meltingMatpro&);

        //- Disallow default bitwise assignment
        void operator=(const UO2meltingMatpro&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2Matpro");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        UO2meltingMatpro
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~UO2meltingMatpro();


    // Member Functions
    
        //- Return true if material is failed
        virtual bool isFailed(const labelList& addr);

        //- Return name of failure criterion
        virtual word criterionName()
        {
            return typeName_();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
