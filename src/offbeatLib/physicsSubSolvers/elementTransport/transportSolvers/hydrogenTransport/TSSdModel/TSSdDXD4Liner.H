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
    Foam::TSSdModel

Description
    TSSD model for the DXD4 outer liner
    From Gong et al.(2019) https://doi.org/10.1016/j.jnucmat.2019.151757

Usage

    \verbatim
    elementTransportOptions
    {
        ...

        HSolverOptions
        {
            ...
            materials
            {
                cladding
                {
                    TSSdModel   TSSdKammenzind;
                    TSSpModel   TSSpKammenzind;
                }
                liner
                {
                    TSSdModel   TSSdDXD4Liner;
                    TSSpModel   TSSpDXD4Liner;
                }
            }
        }
    }
    \endverbatim

SourceFiles
    TSSdDXD4Liner.C

\mainauthor
    Lokesh Verma - Paul Scherrer Institut (PSI)

\contribution

\date 
    July 2024

\*---------------------------------------------------------------------------*/

#ifndef TSSdDXD4Liner_H
#define TSSdDXD4Liner_H

#include "TSSdModel.H"
#include "fvMesh.H"
#include "Table.H"
#include "zeroGradientFvPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "TSSdDXD4Liner" Declaration
\*---------------------------------------------------------------------------*/

class TSSdDXD4Liner
:
    public TSSdModel
{
        //- Disallow default bitwise copy construct
        TSSdDXD4Liner(const TSSdDXD4Liner&);

        //- Disallow default bitwise assignment
        void operator=(const TSSdDXD4Liner&);
        
protected:
        
        // Protected data

public:
 
    //- Runtime type information
    TypeName("TSSdDXD4Liner");

    // Constructors

        //- Construct from dictionary
        TSSdDXD4Liner
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~TSSdDXD4Liner();

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
