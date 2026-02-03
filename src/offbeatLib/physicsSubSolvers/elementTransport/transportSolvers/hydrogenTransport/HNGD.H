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
    Foam::PuRedistributionSolver

Description:
    Hydrides Nucleation Growth Dissolution Model.

    Model the diffusion of hydrogen in solid solution (Css field) due 
    to both the gradient of concentration and of temperature.

    Model the nucleation, growth and dissolution of hydrogen into solid precipitate hydrides.
    Cpp field.

    Parameters of the model (diffusion coefficient, solubility limit etc.) can be tweaked by modifying the C file.
    Current values are took from the refs given belows.

    User can enter:
        Initial concentration in solid solution (flat profile)
        Initial concentration in hydrides form (flat profile)
        Radial profile can be used by modifying the C file. Exemple of such profiles is commented on lines 62-80.
    
        The two following parameters increases the solubility limit as 
        hydrogen precipitates into solid hydrides. Assumes that the precipitation of hydrides deforms
        the matrix so it can accommodate more hydrogen in solution.

        g parameter (default value is 200 wt.ppm). plays a role a low hydride content.
        delta parameter (default value is 1.05). plays a role at high hydride content (no impact below 2000 wt.ppm)

        tau parameter. (default value is 10e4 s). This parameter introduce a time dependance on the precipitation limits.
        During a temperature plateau, the precipitation limit decrease exponentially at a rate of exp( (t-t0) / tau).
        Make the hypotheses that during a temperature hold, H migrates into matrix inhomogeneities reducing as such the 
        precipitation limit.

Usage
    In solverDict:
    \verbatim
    elementTransportOptions
    {
        elements        (H);
    }
    \endverbatim

SourceFiles
    HNGD.C

\mainauthor
    M.Reymond - EPFL

\contribution
    E.Brunetto, A. Scolaro, C. Fiorina - EPFL

\date 
    October 2022

\*---------------------------------------------------------------------------*/

#ifndef HNGD_H
#define HNGD_H

#include "hydrogenTransportSolver.H"
#include "fuelMaterial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class HNGD Declaration
\*---------------------------------------------------------------------------*/

class HNGD
:
    public hydrogenTransportSolver
{
    // Private data
    
        //- Kinetics parameter for dissolution
        volScalarField KD_;
        
        //- Kinetics parameter for nucleation
        volScalarField KN_;
        
        //- Kinetics parameter for growth
        volScalarField KG_;

        //- Solid solution rate coefficient
        volScalarField Kss_;

        //- Explicit source for rate
        volScalarField S_;

        //- Parameter for hydride change
        volScalarField hydrideChange_;

        //- Temperature plateau
        volScalarField t0_;

        //- Enable stress driven diffusion
        bool stressDrivenDiffusion_;
        
        //- Partial molar volume of hydrogen in Zr (m3/mol)
        dimensionedScalar V_H_;
        
    // Private Member Functions

        //- Compute hydride formation energy
        virtual scalar hydrideFormationEnergy(const scalar&);

        //- Convert concentration (wt.ppm) to atomic fraction
        virtual scalar convertToAtomFrac(const scalar&);

        //- compute volume fraction of alpha phase
        virtual scalar fAlpha(const scalar&, const scalar&, const scalar&, const scalar&);

        //- lever rule function
        virtual scalar leverRule(const scalar&, const scalar&, const scalar&);

        //- compute kinetic parameter for growth KG
        virtual scalar KG(const scalar&, const scalar&, const scalar&, 
                          const scalar&);

        //- Update coefficients
        virtual void updateCoefficients();

        //- Update the rates of dissolution, growth and nucleation
        virtual void updateRates();

        //- Disallow default bitwise copy construct
        HNGD(const HNGD&);

        //- Disallow default bitwise assignment
        void operator=(const HNGD&);
     

public:

    //- Runtime type information
    TypeName("HNGD");

    // Constructors

        //- Construct from mesh, materials and dictionary
        HNGD
        (
            fvMesh& mesh,
            const materials& mat,
            const dictionary& elementTransportDict,
            const word& solverName
        );
    
    //- Destructor
    ~HNGD()
    {}

    // Access Functions

    // Member Functions
        
        //- Correct
        virtual void correct();

        virtual void updateMesh()
        {};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
