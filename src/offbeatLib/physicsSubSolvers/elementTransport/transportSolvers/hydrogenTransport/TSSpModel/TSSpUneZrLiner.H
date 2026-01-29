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
    Foam::TSSpModel

Description
    TSSP model for the pure Zirconium inner liner
    From Une et al.(2004) https://doi.org/10.1080/18811248.2004.9715569

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
                    TSSpModel   TSSpUneZr2;
                }
                liner
                {
                    TSSpModel   TSSpUneZrLiner;
                }
            }
        }
    }
    \endverbatim

SourceFiles
    TSSpUneZrLiner.C

\mainauthor
   Lokesh Verma - Paul Scherrer Institut (PSI)

\contribution

\date 
   July 2024

\*---------------------------------------------------------------------------*/

#ifndef TSSpUneZrLiner_H
#define TSSpUneZrLiner_H

#include "TSSpModel.H"
#include "fvMesh.H"
#include "Table.H"
#include "zeroGradientFvPatchField.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class "TSSpUneZrLiner" Declaration
\*---------------------------------------------------------------------------*/

class TSSpUneZrLiner
:
    public TSSpModel
{
        //- Disallow default bitwise copy construct
        TSSpUneZrLiner(const TSSpUneZrLiner&);

        //- Disallow default bitwise assignment
        void operator=(const TSSpUneZrLiner&);
        
protected:
        
        // Protected data

    
public:
 
    //- Runtime type information
    TypeName("TSSpUneZrLiner");

    // Constructors

        //- Construct from dictionary
        TSSpUneZrLiner
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& lawDict
        );

    //- Destructor
   virtual ~TSSpUneZrLiner();

    // Member Functions    
    
    //- Update TSSp
    virtual void correctTSSp
    (
        const labelList& addr
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
