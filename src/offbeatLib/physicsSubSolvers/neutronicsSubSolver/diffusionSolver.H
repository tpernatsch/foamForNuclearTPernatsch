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
    In the `diffusionSolver` solver, the neutron flux distribution is 
    obtained from the solution of the one-group neutron diffusion equation 
    using fitted parameters.

    #### Formulation

    The neutron diffusion solver solves the following equation:

    $$
    \begin{aligned} 
    D\nabla^2 \Phi = \Sigma_A\Phi
    \end{aligned}
    $$

    where:

    - \f$D\f$ is the neutron diffusion coefficient
    - \f$\Phi\f$ is the neutron thermal flux 
    - \f$\Sigma_A\f$ is the macroscopic absorption cross-section

    The obtained thermal neutron flux shape is used by the `Lassmann` burnup class
    to compute the intra-pin power distribution. The solver considers the fuel pin 
    to be a strongly absorbing medium for thermal neutrons, and therefore only diffusion 
    and absorption terms appear in the equation. The boundary condition for neutron 
    thermal flux needed by the solver can be selected by the user to a value of choice. 
    Since only the thermal neutron shape is of importance for the `Lassmann` burnup class,
    the user could simply select the boundary condition for flux to be a `fixedValue` to 1. 
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`diffusion` neutronics solver class:

    \verbatim
        neutronicsSolver diffusion;
    \endverbatim

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
        const dictionary& neutronicsDict
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
