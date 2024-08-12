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
    Foam::ZhangHastelloyCreepModel

Description
    Applies a power law time exponential to compute creep strain rate.
    Model derived from work of Zhang.

SourceFiles
    ZhangHastelloyCreepModel.C

\mainauthor    
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. de Lara - University of Cambridge

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef ZhangHastelloyCreepModel_H
#define ZhangHastelloyCreepModel_H

#include "creepModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class ZhangHastelloyCreepModel
:
    public creepModel
{
    // Private Member Functions

        //- Name used for the fastFluence field
        const word fastFluenceName_;

        //- Name used for the fast flux field
        const word fastFluxName_;
        
        //- Pointer to fastFluence field
        const volScalarField* fastFluence_;
        
        //- Pointer to fast flux field
        const volScalarField* fastFlux_;
        
        //- Index to keep track of time
        scalar timeIndex_;

        //- Relaxation factor
        scalar relax_;

        // Activate Newton-Raphson method for convergence
        Switch NewtonRaphsonMethod_;

        // Tolerance for inner convergence of the Newton-Raphson method
        scalar NewtonRaphsonTolerance_;

        //- Material dependent constant
        scalar A_;

        //- Activation energy
        scalar Qc_;

        //- Multiplier for dpa
        scalar a_;

        //- Exponent for dpa
        scalar b_;

        //- Creep dpa dependence constant
        scalar n0_;

        //- Universal gas constant
        scalar R_;
        
        scalar maxAverageCreep_;
        scalar maxMaximumCreep_;
        
        //- Disallow default bitwise copy construct
        ZhangHastelloyCreepModel(const ZhangHastelloyCreepModel&);

        //- Disallow default bitwise assignment
        void operator=(const ZhangHastelloyCreepModel&);

public:
 
    //- Runtime type information
    TypeName("ZhangHastelloyCreepModel");

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        ZhangHastelloyCreepModel
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    //- Destructor
    ~ZhangHastelloyCreepModel();


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
