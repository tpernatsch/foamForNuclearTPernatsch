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
    Foam::diffusionSolver

Description
    Daugther class of neutronicsSubSolver class based on one-group diffusion 
    theory. 

    \htmlinclude neutronicsSubSolver_diffusionSolver.html

SourceFiles
    diffusionSolver.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n

\contribution
    A. Scolaro, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef diffusionSolver_H
#define diffusionSolver_H

#include "neutronicsSubSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class diffusionSolver Declaration
\*---------------------------------------------------------------------------*/

class diffusionSolver
:
public neutronicsSubSolver
{
    // Private data

    // Private Member Functions

protected:

    //- Thermal neutron flux distribution
    volScalarField neutronFlux0_;

    //- Square inverse diffusion length
    volScalarField invDiffLength_;

    //- Macroscopic absorption cross section
    volScalarField SigmaA_;

    //- Neutron diffusion coefficient
    volScalarField neutronD_;

    //- Disallow default bitwise copy construct
    diffusionSolver(const diffusionSolver&);

    //- Disallow default bitwise assignment
    void operator=(const diffusionSolver&);
     

public:

    //- Runtime type information
    TypeName("diffusion");

    // Constructors

    //- Construct from mesh, bu and dictionary
    diffusionSolver
    (
        const fvMesh& mesh,
        const burnup& bu,
        const dictionary& neutronicsOptDict
    );
    
    //- Destructor
    ~diffusionSolver()
    {}

    // Access functions

    //- Return a const reference to the neutron flux field 
    inline const volScalarField& neutronFlux0() const
    {
        return neutronFlux0_;
    }

    // Member Functions
        
    //- Update neutronFlux0
    virtual void correct();

    //- Return true if converged
    virtual bool converged();    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
