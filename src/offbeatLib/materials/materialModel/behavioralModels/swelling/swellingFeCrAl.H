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
    Foam::swellingFeCrAl

Description
    Class handling swelling phenomenon for FeCrAl derived from Bison manual.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            swellingModel     FeCrAl;
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingFeCrAl.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef swellingFeCrAl_H
#define swellingFeCrAl_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingFeCrAl Declaration
\*---------------------------------------------------------------------------*/

class swellingFeCrAl
:
    public swellingModel
{
    // Private data

        //- Name used for the burnup field
        const word fastFluenceName_;
        
        //- Reference to Burnup field
        const volScalarField* fastFluence_;

        //- Parameters for the FeCrAl model
        scalar par1;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingFeCrAl(const swellingFeCrAl&);

        //- Disallow default bitwise assignment
        void operator=(const swellingFeCrAl&);
    
    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(swellingFeCrAl, group_, group);
    
        //- Swelling strain scaling factor [-]
        declareParameterWithVirtual(swellingFeCrAl, 
            F_epsilonSwelling_, F_epsilonSwelling);
        
        //- Swelling strain offset [-]
        declareParameterWithVirtual(swellingFeCrAl, 
            delta_epsilonSwelling_, delta_epsilonSwelling);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("FeCrAl");

    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        swellingFeCrAl
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~swellingFeCrAl();


    // Member Functions
    
        //- Update heatCapacity  
        virtual void correct(const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
