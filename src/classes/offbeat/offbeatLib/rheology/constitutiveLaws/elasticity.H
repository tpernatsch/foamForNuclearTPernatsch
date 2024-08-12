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
    Foam::elasticity

Description
    Daugther class of constitutveLaw class that implements the liner Hooke's law mechanical behavior.

    \htmlinclude constitutiveLaw_elasticity.html

SourceFiles
    elasticity.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef elasticity_H
#define elasticity_H

#include "constitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermomechanicsSubSolver Declaration
\*---------------------------------------------------------------------------*/

class elasticity
:
    public constitutiveLaw
{
    // Private data  


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        elasticity(const elasticity&);

        //- Disallow default bitwise assignment
        void operator=(const elasticity&);


public:
 
    //- Runtime type information
    TypeName("elasticity");


    // Constructors

        //- Construct from mesh, materials , dict and labelList
        elasticity
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    ~elasticity();


    // Member Functions

        //- Update mechanical law     
        virtual void correct
        (
            volScalarField& sigmaHyd, 
            volSymmTensorField& sigmaDev, 
            volSymmTensorField& epsEl,
            const labelList& addr
        );
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
