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
    Foam::swellingGrowthMatproZy

Description
    Class modelling the irradiation growth phenomenon according to MATPRO 
    correlations.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            swellingModel     growthMATPROZy;
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingGrowthMatproZy.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef swellingGrowthMatproZy_H
#define swellingGrowthMatproZy_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingGrowthMatproZy Declaration
\*---------------------------------------------------------------------------*/

class swellingGrowthMatproZy
:
    public swellingModel
{
    // Private data

        //- Name used for the fastFluence field
        const word fastFluenceName_;
        
        //- Reference to fastFluence field
        const volScalarField& fastFluence_;

        //- Cold work
        scalar cW_;

        //- Parameters for the Matpro model
        scalar A;
        scalar par1;
        scalar fz;

        //- Perturbation parameter
        scalar perturb;

        // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingGrowthMatproZy(const swellingGrowthMatproZy&);

        //- Disallow default bitwise assignment
        void operator=(const swellingGrowthMatproZy&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("growthMATPROZy");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        swellingGrowthMatproZy
        (
            const fvMesh& mesh,
            const dictionary& dict 
        );


    // Selectors


    //- Destructor
    virtual ~swellingGrowthMatproZy();


    // Member Functions
    
        //- Update irradiation growth
        virtual void correct(const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
