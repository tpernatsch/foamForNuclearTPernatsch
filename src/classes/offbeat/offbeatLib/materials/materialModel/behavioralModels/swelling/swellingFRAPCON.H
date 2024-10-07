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
    Foam::swellingFRAPCON

Description
    Class modelling swelling phenomenon derived from FRAPCON.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            swellingModel     UO2FRAPCON;
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingFRAPCON.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef swellingFRAPCON_H
#define swellingFRAPCON_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingFRAPCON Declaration
\*---------------------------------------------------------------------------*/

class swellingFRAPCON
:
    public swellingModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;

        //- Name used for the density field
        const word densityName_;
        
        //- Reference to density field
        const volScalarField* rho_;

        //- Parameters for the FRAPCON model
        scalar par1;
        scalar par2;
        scalar par3;
        scalar par4;
        scalar par5;
        scalar par6;

        //- Perturbation parameter
        scalar perturb;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingFRAPCON(const swellingFRAPCON&);

        //- Disallow default bitwise assignment
        void operator=(const swellingFRAPCON&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UO2FRAPCON");

    // Declare run-time constructor selection table

    // Constructors

        //- Construct from dictionary
        swellingFRAPCON
        (
            const fvMesh& mesh,
            const dictionary& dict  
        );


    // Selectors

    //- Destructor
    virtual ~swellingFRAPCON();


    // Member Functions
    
    //- Update heatCapacity  
    virtual void correct(const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
