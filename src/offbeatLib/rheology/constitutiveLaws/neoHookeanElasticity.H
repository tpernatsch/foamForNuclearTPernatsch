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
    Foam::neoHookeanElasticity

Description
    NeoHookean elasticity constitutive law.
    
    \htmlinclude constitutiveLaw_neoHookeanElasticity.html

SourceFiles
    neoHookeanElasticity.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef neoHookeanElasticity_H
#define neoHookeanElasticity_H

#include "constitutiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class thermomechanicsSubSolver Declaration
\*---------------------------------------------------------------------------*/

class neoHookeanElasticity
:
    public constitutiveLaw
{
    // Private data  


    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        neoHookeanElasticity(const neoHookeanElasticity&);

        //- Disallow default bitwise assignment
        void operator=(const neoHookeanElasticity&);


public:
 
    //- Runtime type information
    TypeName("neoHookeanElasticity");


    // Constructors

        //- Construct from mesh, materials , dict and labelList
        neoHookeanElasticity
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );


    //- Destructor
    ~neoHookeanElasticity();


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
