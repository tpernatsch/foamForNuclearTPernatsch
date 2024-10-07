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
    Foam::powerLaw

Description
    Applies a power law to compute creep strain rate.

    \htmlinclude creepModel_powerLaw.html

SourceFiles
    powerLaw.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef powerLaw_H
#define powerLaw_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class powerLaw
:
    public creepModel
{
    // Private Member Functions
        
        scalar timeIndex_;

        scalar relax_;
        
        scalar maxAverageCreep_;
        scalar maxMaximumCreep_;

        //- Dimensionless constant
        const scalar B; 
        const scalar sigmaC; 
        const scalar n;
        
        //- Disallow default bitwise copy construct
        powerLaw(const powerLaw&);

        //- Disallow default bitwise assignment
        void operator=(const powerLaw&);

public:
 
    //- Runtime type information
    TypeName("powerLaw");

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        powerLaw
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~powerLaw();


    // Member Functions

    
        //- Update mechanical properties according to rheology model        
        virtual void correctCreep
        (
            volSymmTensorField& epsilonEl,
            const labelList& addr
        );
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
