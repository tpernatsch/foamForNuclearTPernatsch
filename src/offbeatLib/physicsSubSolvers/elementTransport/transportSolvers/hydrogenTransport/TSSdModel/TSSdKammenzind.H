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
    Foam::creepModel

Description
    varying yield stress as a function of plastic strain

Usage

    \verbatim
    rheology byMaterial;
    
    ...

    // List of materials, one per cellZone.
    materials
    {
        ...

        clad
        {
            material Zircaloy;
            ...

            rheologyModel misesPlasticity;

            rheologyModelOptions
            {   
                yieldStressModel    TSSdKammenzind;
                plasticStrainVsYieldStress
                (
                    (0.0 200e6)
                    (0.005 250e6)
                    (0.01 300e6)
                );
                creepModel	none;
            }
        }
    }
    \endverbatim

SourceFiles
    yieldStress.C

\mainauthor
    Matthieu Reymond - EPFL

\contribution

\date 
    Feb. 2023

\*---------------------------------------------------------------------------*/

#ifndef TSSdKammenzind_H
#define TSSdKammenzind_H

#include "TSSdModel.H"
#include "fvMesh.H"
#include "Table.H"
#include "zeroGradientFvPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "TSSdKammenzind" Declaration
\*---------------------------------------------------------------------------*/

class TSSdKammenzind
:
    public TSSdModel
{
        //- Disallow default bitwise copy construct
        TSSdKammenzind(const TSSdKammenzind&);

        //- Disallow default bitwise assignment
        void operator=(const TSSdKammenzind&);
        
protected:
        
        // Protected data
public:
 
    //- Runtime type information
    TypeName("TSSdKammenzind");

    // Constructors

        //- Construct from dictionary
        TSSdKammenzind
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~TSSdKammenzind();

    // Member Functions    
    
    //- Update TSSd
    virtual void correctTSSd
    (
        const labelList& addr
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
