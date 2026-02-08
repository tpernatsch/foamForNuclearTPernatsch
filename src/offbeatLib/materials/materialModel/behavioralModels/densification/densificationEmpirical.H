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
    Foam::densificationEmpirical

Description
    Empirical model for densification. The fuel tends toward a final value of 
    density (due to densification) with an exponential function of Bu.

    The user can introduce the final density increase due to densification 
    and the time constant of the densification process.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UO2;

            ...

            densificationModel         Empirical;
            
            // Final density increase (in %) once densification is complete
            densificationDensityChange   0.5;
            
            // Time constant (in MWd/tU) for densification process
            densificationTimeConstant   1000;
        }

        ...
    }
    \endverbatim

SourceFiles
    densificationEmpirical.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef densificationEmpirical_H
#define densificationEmpirical_H

#include "densificationModel.H"
#include "physicoChemicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class densificationEmpirical Declaration
\*---------------------------------------------------------------------------*/

class densificationEmpirical
:
    public densificationModel
{
    // Private data

        //- Name used for the burnup field
        const word burnupName_;
        
        //- Reference to Burnup field
        const volScalarField* Bu_;
        
        //- Fuel density as fraction of TD
        scalar densityFrac_;        
        
        //- Final density change once densification is complete
        scalar densificationDensityChange_;

        //- Time constant (in MWd/tU) determining the speed of densification
        scalar densificationTimeConstant_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        densificationEmpirical(const densificationEmpirical&);

        //- Disallow default bitwise assignment
        void operator=(const densificationEmpirical&);

    // Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(densificationEmpirical, group_, group);
    
        //- Densification strain scaling factor [-]
        declareParameterWithVirtual(densificationEmpirical, 
            F_epsilonDensification_, F_epsilonDensification);
        
        //- Densification strain offset [-]
        declareParameterWithVirtual(densificationEmpirical, 
            delta_epsilonDensification_, delta_epsilonDensification);


protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("empirical");

    
    // Declare run-time constructor selection table

    
    // Constructors

        //- Construct from dictionary
        densificationEmpirical
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );

    // Selectors

        //- Destructor
        virtual ~densificationEmpirical();


    // Member Functions
    
        //- Update heatCapacity  
        virtual void correct(const scalarField& T, const labelList& addr);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
