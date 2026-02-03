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
    Foam::constantSwelling

Description
    Class caculating the swelling with a provided constant swelling rate.

    \htmlinclude swelling_constantSwelling.html

SourceFiles
    constantSwelling.C

\mainauthor
    F. Xiang - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour) & XJTU (Xi'an JiaoTong
    University, China, Nuclear Thermal-hydraulic Laboratory)

\contribution
    A. Scolaro - EPFL
    E. Brunetto - EPFL

\date
    January 2023

\*---------------------------------------------------------------------------*/

#ifndef constantSwelling_H
#define constantSwelling_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantSwelling Declaration
\*---------------------------------------------------------------------------*/

class constantSwelling
:
    public swellingModel
{
    // Private data

        //- Swelling rate. In unit of /(10^25n/m²)
        scalar s_;

        //- Name used for the burnup field
        const word fastFluenceName_;

        //- Reference to fast neutron fluence field
        const volScalarField* fastFluence_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantSwelling(const constantSwelling&);

        //- Disallow default bitwise assignment
        void operator=(const constantSwelling&);


    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(constantSwelling, group_, group);
    
        //- Swelling strain scaling factor [-]
        declareParameterWithVirtual(constantSwelling, 
            F_epsilonSwelling_, F_epsilonSwelling);
        
        //- Swelling strain offset [-]
        declareParameterWithVirtual(constantSwelling, 
            delta_epsilonSwelling_, delta_epsilonSwelling);

protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("constantRate");

    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        constantSwelling
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~constantSwelling();


    // Member Functions

        //- Update
        virtual void correct(const scalarField& T, const labelList& addr) ;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
