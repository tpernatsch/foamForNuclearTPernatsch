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
    Foam::poreVelocityUO2Sens

Description
    Model to compute pore velocity with the correlation provided by Sens (P.F.
    Sens. The kinetics of pore movement in UO2 fuel rods. Journal of Nuclear 
    Materials, 43(3):293–307, 1972. doi:10.1016/0022-3115(72)90061-X.).

    In addition, a correction factor can be defined to correct the velocity 
    (default value is 1.0).

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UPuO2;

            ...

            poreVelocityModel     UO2Sens;
            // Multiplicator for pore velocity
            // alphaPoreVelocity 1;
        }

        ...
    }
    \endverbatim

SourceFiles
    poreVelocityUO2Sens.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef poreVelocityUO2Sens_H
#define poreVelocityUO2Sens_H

#include "poreVelocityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poreVelocityUO2Sens Declaration
\*---------------------------------------------------------------------------*/

class poreVelocityUO2Sens
:
    public poreVelocityModel
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        poreVelocityUO2Sens(const poreVelocityUO2Sens&);

        //- Disallow default bitwise assignment
        void operator=(const poreVelocityUO2Sens&);

protected:

        //- Pointer to temperature field
        const volScalarField* T_;

        //- Pointer to temperature gradient field
        const volVectorField* gradT_;

        //- Universal gas constant (erg/mol/K)
        const scalar R_;

public:

    //- Runtime type information
    TypeName("UO2Sens");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        poreVelocityUO2Sens
        (
            const fvMesh& mesh, 
            const dictionary& dict,
            const word defaultModel
        );


    //- Destructor
    ~poreVelocityUO2Sens();

    // Access functions    

    // Member Functions

        //- Correct velocity according to Sens model
        virtual void correct( const labelList& addr );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
