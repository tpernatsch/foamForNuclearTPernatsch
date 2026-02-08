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
    Foam::HTransportSolver

Description:
    This solver calculate the diffusion of hydrogen in solid solution in the cladding
    due to a temperature or concentration gradient.
    It also calculates the precipitation or dissolution of the hydrogen into solid hydrides.

    3 fields are considered: - Concentration in solid solution Css (wt.ppm)
                                      - Concentration in precipitates Cpp (wt.ppm)
                                      - Total concentration of H Ctot (wt.ppm)

    The user must provide the initial Css and Cpp fields.
    For Css field, use the hydrogenFvPatchField BC (see documentation)

    This solvers considers the supersolubility limit (TSSp) and the solubiltiy limit (TSSd)
    of hydrogen in zyrcalloy.
    4 cases are thus possible:
                                - precipitation
                                - dissolution
                                - phase equilibrium

    The kinetic of precipitation is considered by introducing the kinetic parameter
    as proposed by Kammenzind et al. It can be switched off for the simulation of diffusion driven problems
    during which the kinetic of precipitation is negligeable.
    

Usage
    In solverDict:
    \verbatim
    elementTransportOptions
    {
        solvers
        (
            HRedistribution
        );
        HSolverOptions
        {
            precipitationKinetic off; //- switch off the kinetic of precipitation for diffusion driven problems
                                      //- set to on by default
            hardCodedInitialProfile on;  //- advanced users. possibility to hard code the initial profiles (Css/Cpp). Set to off by default.
                                      //- in the solver. The use of codeStream in the /0/Css and /0/Cpp is recommanded for custom initial profiles
        }
    }

    \endverbatim

SourceFiles
    HSolver.C

\mainauthor
    M.Reymond - EPFL

\contribution
    E.Brunetto, A. Scolaro, C. Fiorina - EPFL
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef HSolver_H
#define HSolver_H

#include "hydrogenTransportSolver.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class HSolver Declaration
\*---------------------------------------------------------------------------*/

class HSolver
:
public hydrogenTransportSolver
{
    // Private data
        //- kinetic parameter of precipitation [s-0.5]
        volScalarField alpha_;
        
        scalar kineticParameter_;

    // Private Member Functions

        //- Update diffusion coefficient and TSS
        virtual void updateCoefficients();
        
        //- Update the rate terms of precipitation/dissolution
        virtual void updateRateTerms();

        
        //- Disallow default bitwise copy construct
        HSolver(const HSolver&);

        //- Disallow default bitwise assignment
        void operator=(const HSolver&);
     

public:

    //- Runtime type information
    TypeName("HSolver");

    // Constructors

        //- Construct from mesh, materials and dictionary
        HSolver
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict
        );
    
    //- Destructor
    ~HSolver()
    {}

    // Access Functions

    // Return name of solving variable
    virtual const word& fieldName()
    {
        return Css_.name();
        return Cpp_.name();
    };

    // Member Functions
        
    //- Correct
    virtual void correct();

    virtual void updateMesh(){};
    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
