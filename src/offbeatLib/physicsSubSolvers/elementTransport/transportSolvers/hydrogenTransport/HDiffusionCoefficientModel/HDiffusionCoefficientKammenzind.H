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
    Foam::HDiffusionCoefficientKammenzind

Description
    Diffusion coefficient of Hydrogen in Zircaloy-4 from Kammenzind et al.
    
    W.J. Kammenzind, B. Franklin, D.G. Peters, H.R. Dufﬁn, 
    Hydrogen pickup and redistribution in alpha-annealed Zircaloy-4,
    11th International Symposium on Zirconium in the Nuclear Industry, 
    ASTM STP1295 (1997), p. 338, 10.2172/10181237

Usage

SourceFiles
    HDiffusionCoefficientKammenzind.C

\mainauthor
    Matthieu Reymond - EPFL

\contribution

\date 
    Feb. 2023

\*---------------------------------------------------------------------------*/

#ifndef HDiffusionCoefficientKammenzind_H
#define HDiffusionCoefficientKammenzind_H

#include "HDiffusionCoefficientModel.H"
#include "fvMesh.H"
#include "Table.H"
#include "zeroGradientFvPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "HDiffusionCoefficientKammenzind" Declaration
\*---------------------------------------------------------------------------*/

class HDiffusionCoefficientKammenzind
:
    public HDiffusionCoefficientModel
{
        //- Disallow default bitwise copy construct
        HDiffusionCoefficientKammenzind(const HDiffusionCoefficientKammenzind&);

        //- Disallow default bitwise assignment
        void operator=(const HDiffusionCoefficientKammenzind&);
        
protected:
        
        // Protected data

public:
 
    //- Runtime type information
    TypeName("HDiffusionCoefficientKammenzind");

    // Constructors

        //- Construct from dictionary
        HDiffusionCoefficientKammenzind
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~HDiffusionCoefficientKammenzind();

    // Member Functions    
    
    //- Update TSSd
    virtual void correctDH
    (
        const labelList& addr
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
