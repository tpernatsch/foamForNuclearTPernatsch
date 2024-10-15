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
    Foam::swellingWrightShamHastelloyN

Description
    Class modelling the void swelling phenomenon for HastelloyN cladding
    material according to Wright and Sham correlations.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material HastelloyN;

            ...

            // WrightSham is already selected by default for HastelloyN material
            swellingModel     WrightShamHastelloyN;
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingWrightShamHastelloyN.C

 \mainauthor
    Alejandra de Lara - University of Cambridge

 \contribution
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)\n

 \date
    Sep 2023

\*---------------------------------------------------------------------------*/

#ifndef swellingWrightShamHastelloyN_H
#define swellingWrightShamHastelloyN_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingWrightShamHastelloyN Declaration
\*---------------------------------------------------------------------------*/

class swellingWrightShamHastelloyN
:
    public swellingModel
{
    // Private data

        //- Name used for the fastFluence field
        const word fastFluenceName_;
        
        //- Reference to fastFluence field
        const volScalarField* fastFluence_;

        //- Parameters for the generalized model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingWrightShamHastelloyN(const swellingWrightShamHastelloyN&);

        //- Disallow default bitwise assignment
        void operator=(const swellingWrightShamHastelloyN&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("WrightShamHastelloyN");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        swellingWrightShamHastelloyN
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors


    //- Destructor
    virtual ~swellingWrightShamHastelloyN();


    // Member Functions
    
        //- Update irradiation growth
        virtual void correct(const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

