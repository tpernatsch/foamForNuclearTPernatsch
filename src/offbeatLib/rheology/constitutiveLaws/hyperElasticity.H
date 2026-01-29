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
    Foam::hyperElasticity

Description
    StVenant-Kirchhoff hyperElasticity mechanical behaviour.
    
    \htmlinclude constitutiveLaw_hyperElasticity.html

SourceFiles
    hyperElasticity.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef hyperElasticity_H
#define hyperElasticity_H

#include "constitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermomechanicsSubSolver Declaration
\*---------------------------------------------------------------------------*/

class hyperElasticity
:
    public constitutiveLaw
{
    // Private data  


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        hyperElasticity(const hyperElasticity&);

        //- Disallow default bitwise assignment
        void operator=(const hyperElasticity&);


public:
 
    //- Runtime type information
    TypeName("hyperElasticity");


    // Constructors

        //- Construct from mesh, materials , dict and labelList
        hyperElasticity
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~hyperElasticity();


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
