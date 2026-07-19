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
    Foam::poreVelocityMOXClementFinnis

Description
    Model to compute pore velocity with the correlation provided by Clement and
    Finnis.

    Ref: 
    C.F. Clement, M.W. Finnis,
    Plutonium redistribution in mixed oxide (U, Pu)O2 nuclear fuel elements,
    Journal of Nuclear Materials,
    Volume 75, Issue 1,
    1978,
    Pages 193-200,
    ISSN 0022-3115,
    https://doi.org/10.1016/0022-3115(78)90044-2.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        fuel
        {
            material UPuO2;

            ...

            poreVelocityModel     MOXClementFinnis;
            // Multiplicator for pore velocity
            // alphaPoreVelocity 1;
        }

        ...
    }
    \endverbatim

    In addition, a correction factor can be defined to correct the velocity 
    (default value is 1.0).

SourceFiles
    poreVelocityMOXClementFinnis.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef poreVelocityMOXClementFinnis_H
#define poreVelocityMOXClementFinnis_H

#include "poreVelocityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poreVelocityMOXClementFinnis Declaration
\*---------------------------------------------------------------------------*/

class poreVelocityMOXClementFinnis
:
    public poreVelocityModel
{
    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        poreVelocityMOXClementFinnis(const poreVelocityMOXClementFinnis&);

        //- Disallow default bitwise assignment
        void operator=(const poreVelocityMOXClementFinnis&);

protected:

        //- Pointer to temperature field
        const volScalarField* T_;

        //- Pointer to temperature gradient field
        const volVectorField* gradT_;

        //- Constant for correlation (m/s)
        const scalar A_;

        //- Constant for correlation (K)
        const scalar alpha_;

        //- Pore diameter
        const scalar l_;
 

public:

    //- Runtime type information
    TypeName("MOXClementFinnis");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        poreVelocityMOXClementFinnis
        (
            const fvMesh& mesh, 
            const dictionary& dict,
            const word defaultModel
        );


    //- Destructor
    ~poreVelocityMOXClementFinnis();

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
