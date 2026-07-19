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
    Foam::poreVelocityMOXVaporPressure

Description

Usage


SourceFiles
    poreVelocityMOXVaporPressure.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2022

\*---------------------------------------------------------------------------*/

#ifndef poreVelocityMOXVaporPressure_H
#define poreVelocityMOXVaporPressure_H

#include "poreVelocityModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class poreVelocityMOXVaporPressure Declaration
\*---------------------------------------------------------------------------*/

class poreVelocityMOXVaporPressure
:
    public poreVelocityModel
{
    // Private data

        //- Pointer to temperature field
        const volScalarField* T_;

        //- Pointer to temperature gradient field
        const volVectorField* gradT_;

        //- Pointer to oxygen to metal ratio field
        const volScalarField* OM_;

        //- PuO2 lattice parameter (cm)
        const scalar CPuO2_;

        //- UO2 lattice parameter (cm)
        const scalar CUO2_;

        //- Pu mass fraction
        scalar wPu_;

        //- U mass fraction
        scalar wU_;

        //- Am mass fraction
        scalar wAm_;

        //- PuO2 molar fraction
        scalar xPuO2_;

        //- UO2 molar fraction
        scalar xUO2_;

        //- Boltzmann constant (cm2 * g / s / K)
        const scalar kb_;

        //- Gas molecules' diameter (cm)
        const scalar d_;

        //- Sum of molecular masses of N2 and UO2
        const scalar m_;

        //- n constant for diffusion coefficient 
        // (value from fig.A.1 of https://doi.org/10.1080/00295450.2017.1314748)
        const scalar nD_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        poreVelocityMOXVaporPressure(const poreVelocityMOXVaporPressure&);

        //- Disallow default bitwise assignment
        void operator=(const poreVelocityMOXVaporPressure&);

        //- Calculate partial pressure oxygen
        virtual const Foam::scalar calcPartialPressureO2(const scalar OMi, const scalar Ti);

        //- Calculate integral for Gibbs' energy 
        virtual const Foam::scalar calcIntegral( const scalar upperLimit, const scalar Ti, const scalar q);

        //- Calculate total pressure
        virtual Foam::scalar calcVaporPressure(const scalar OMi, const scalar Ti);

        //- Calculate total pressure derivative
        virtual Foam::scalar calcVaporPressureDerivative(const scalar OMi, const scalar Ti);


public:

    //- Runtime type information
    TypeName("vaporPressureMOX");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        poreVelocityMOXVaporPressure
        (
            const fvMesh& mesh, 
            const dictionary& dict,
            const word defaultModel
        );


    //- Destructor
    ~poreVelocityMOXVaporPressure();

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
