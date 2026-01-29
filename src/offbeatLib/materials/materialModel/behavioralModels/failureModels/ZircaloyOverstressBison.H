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
    Foam::ZircaloyOverstressBison

Description
    Burst overstress criterion for Zircaloy-4 derived from Bison code. NOTE: 
    oxidation contribution still needs to be included.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        cladding
        {
            failureModel    ZircaloyOverstressBison;
            stopIfFailed    on;
            failureModelOptions
            {
                patchNames      ( "cladOuter" "cladInner" );
            }
        }
    }
    \endverbatim

SourceFiles
    ZircaloyOverstressBison.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    A. Scolaro, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef ZircaloyOverstressBison_H
#define ZircaloyOverstressBison_H

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class ZircaloyOverstressBison Declaration
\*---------------------------------------------------------------------------*/

class ZircaloyOverstressBison
:
    public failureModel
{
    // Private data

        //- List of patch names where the criterion applies
        const wordList patchNames_;

        //- Reference to betaFraction field
        const volScalarField* betaFrac_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        ZircaloyOverstressBison(const ZircaloyOverstressBison&);

        //- Disallow default bitwise assignment
        void operator=(const ZircaloyOverstressBison&);

protected:
    
    // Protected data

public:

    //- Runtime type information
    TypeName("ZircaloyOverstressBison");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from dictionary
        ZircaloyOverstressBison
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const dictionary& baseDict
        );


    // Selectors


    //- Destructor
    virtual ~ZircaloyOverstressBison();


    // Member Functions
    
        //- Return true if material is failed
        virtual bool isFailed(const labelList& addr);

        //- Return name of failure criterion
        virtual word criterionName()
        {
            return typeName_();
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
