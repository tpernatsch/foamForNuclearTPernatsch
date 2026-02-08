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
    Foam::poreVelocityMOXLackey

Description
    Model to compute pore velocity within MOX fuel with the correlation provided
    by Lackey et Al.

    Ref: 
    W. J. Lackey, F. J. Homan & A. R. Olsen (1972) Porosity and Actinide 
    Redistribution During Irradiation of (U, Pu)O2, Nuclear Technology,
    16:1, 120-142, DOI: 10.13182/NT72-A31181

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UPuO2;

            ...

            poreVelocityModel     MOXLackey;
            // Multiplicator for pore velocity
            // alphaPoreVelocity 1;
        }

        ...
    }
    \endverbatim

SourceFiles
    poreVelocityMOXLackey.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef poreVelocityMOXLackey_H
#define poreVelocityMOXLackey_H

#include "poreVelocityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poreVelocityMOXLackey Declaration
\*---------------------------------------------------------------------------*/

class poreVelocityMOXLackey
:
    public poreVelocityModel
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        poreVelocityMOXLackey(const poreVelocityMOXLackey&);

        //- Disallow default bitwise assignment
        void operator=(const poreVelocityMOXLackey&);

protected:

        //- Pointer to temperature field
        const volScalarField* T_;

        //- Pointer to temperature gradient field
        const volVectorField* gradT_;

        //- Pointer to oxygen to metal ratio field
        const volScalarField* OM_;

        //- Correlation factor
        const scalar K_;

public:

    //- Runtime type information
    TypeName("MOXLackey");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        poreVelocityMOXLackey
        (
            const fvMesh& mesh, 
            const dictionary& dict,
            const word defaultModel
        );


    //- Destructor
    ~poreVelocityMOXLackey();

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
