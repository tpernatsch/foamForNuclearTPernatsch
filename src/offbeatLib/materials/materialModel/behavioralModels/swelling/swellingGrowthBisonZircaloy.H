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
    Foam::swellingGrowthBisonZircaloy

Description
    Class modelling the irradiation growth phenomenon derived from Moose documentation.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            swellingModel     ZircaloyBison;

            swellingOptions
            {
                cladType SRA;//RXA;//PRA;//ZIRLO;//ESCORE;//M5;//
            }
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingGrowthBisonZircaloy.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef swellingGrowthBisonZircaloy_H
#define swellingGrowthBisonZircaloy_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingGrowthBisonZircaloy Declaration
\*---------------------------------------------------------------------------*/

class swellingGrowthBisonZircaloy
:
    public swellingModel
{
    // Private data

        //- Name used for the fastFluence field
        const word fastFluenceName_;
        
        //- Reference to fastFluence field
        const volScalarField* fastFluence_;

        //- Parameter for the fast flux model (constant in formula A*phi**n)
        scalar A;

        //- Parameter for the fast flux model (exponent in formula A*phi**n)
        scalar n;

        //- Type of cladding material for growth correlations. Default is ESCORE 
        //  (for consistency with past validation cases. To be changed to SRA)
        word cladType_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingGrowthBisonZircaloy(const swellingGrowthBisonZircaloy&);

        //- Disallow default bitwise assignment
        void operator=(const swellingGrowthBisonZircaloy&);

    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(swellingGrowthBisonZircaloy, group_, group);
    
        //- Swelling strain scaling factor [-]
        declareParameterWithVirtual(swellingGrowthBisonZircaloy, 
            F_epsilonSwelling_, F_epsilonSwelling);
        
        //- Swelling strain offset [-]
        declareParameterWithVirtual(swellingGrowthBisonZircaloy, 
            delta_epsilonSwelling_, delta_epsilonSwelling);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZircaloyBison");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        swellingGrowthBisonZircaloy
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~swellingGrowthBisonZircaloy();


    // Member Functions
    
        //- Update irradiation growth
        virtual void correct(const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
