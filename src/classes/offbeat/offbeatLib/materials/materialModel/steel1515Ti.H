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
    Foam::steel1515Ti

Description
    Class modelling the austenitic steel 15-15 Ti cladding material.

SourceFiles
    steel1515Ti.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef steel1515Ti_H
#define steel1515Ti_H

#include "materialModel.H"
#include "volFields.H"

#include "swellingModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class steel1515Ti Declaration
\*---------------------------------------------------------------------------*/

class steel1515Ti
:
    public materialModel
{
    // Private data

        //- irradiationGrowth model
        autoPtr<swellingModel> swelling_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        steel1515Ti(const steel1515Ti&);

        //- Disallow default bitwise assignment
        void operator=(const steel1515Ti&);

public:

    //- Runtime type information
    TypeName("steel1515Ti");


    // Constructors

        //- Construct from mesh, global options, materials, and dictionary
        steel1515Ti
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~steel1515Ti();


    // Member Functions
    
        //- Update the behavioral models
        virtual void correctBehavioralModels( const scalarField& T );
        
        //- Return the tensor of additional strains
        virtual void additionalStrains( symmTensorField& sf ) const;

        //- Return name of the material
        virtual word name() const
        {
            return "steel1515Ti";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
