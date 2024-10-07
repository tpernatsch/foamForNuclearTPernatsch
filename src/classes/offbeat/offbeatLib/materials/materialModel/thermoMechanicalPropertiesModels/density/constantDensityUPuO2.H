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
    Foam::constantDensityUPuO2

Description
    Class modelling the constant density of (U,Pu)O2 MOX fuel.

SourceFiles
    constantDensityUPuO2.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantDensityUPuO2_H
#define constantDensityUPuO2_H

#include "densityModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantDensityUPuO2 Declaration
\*---------------------------------------------------------------------------*/

class constantDensityUPuO2
:
    public densityModel
{
    // Private data

        //- Fraction of Theoretical Density 
        scalar densityFrac_;

        //- Theoretical density 
        scalar theoreticalDensity_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantDensityUPuO2(const constantDensityUPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const constantDensityUPuO2&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("UPuO2Constant");

    // Constructors

        //- Construct from dictionary
        constantDensityUPuO2
        (
            const fvMesh& mesh,
            const dictionary& dict, 
            const word defaultModel   
        );

    //- Destructor
    virtual ~constantDensityUPuO2();


    // Member Functions
    
    //- Update conductivity  
    virtual void correct(scalarField& sf, const scalarField& T, const labelList& addr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
