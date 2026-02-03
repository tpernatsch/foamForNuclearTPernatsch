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
    Foam::swellingGrowthGeneralized1515Ti

Description
    Class modelling the void swelling growth phenomenon for 15-15 Ti cladding 
    material according to the 'Generalized 15-15 Ti' correlation.
    Source : "Modeling and Analysis of Nuclear Fuel Pin Behavior for Innovative
    Lead Cooled FBR - L.Luzzi et al. - 2014"

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            swellingModel     Steel1515TiGeneralized;
        }

        ...
    }
    \endverbatim

SourceFiles
    swellingGrowthGeneralized1515Ti.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef swellingGrowthGeneralized1515Ti_H
#define swellingGrowthGeneralized1515Ti_H

#include "swellingModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class swellingGrowthGeneralized1515Ti Declaration
\*---------------------------------------------------------------------------*/

class swellingGrowthGeneralized1515Ti
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
        scalar par5;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        swellingGrowthGeneralized1515Ti(const swellingGrowthGeneralized1515Ti&);

        //- Disallow default bitwise assignment
        void operator=(const swellingGrowthGeneralized1515Ti&);

    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(swellingGrowthGeneralized1515Ti, group_, group);
    
        //- Swelling strain scaling factor [-]
        declareParameterWithVirtual(swellingGrowthGeneralized1515Ti, 
            F_epsilonSwelling_, F_epsilonSwelling);
        
        //- Swelling strain offset [-]
        declareParameterWithVirtual(swellingGrowthGeneralized1515Ti, 
            delta_epsilonSwelling_, delta_epsilonSwelling);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("Steel1515TiGeneralized");

   
    // Declare run-time constructor selection table

   
    // Constructors

        //- Construct from dictionary
        swellingGrowthGeneralized1515Ti
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~swellingGrowthGeneralized1515Ti();


    // Member Functions
    
        //- Update irradiation growth
        virtual void correct(const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
