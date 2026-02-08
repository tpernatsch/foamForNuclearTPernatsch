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
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fastFlux

Description
	As the name implies, the `none` fast flux class in OFFBEAT allows you to switch
	off the fast flux and fast fluence models and have a simulation where these
	fields are not considered.
	This means that the fast flux and fast fluence fields `fastFlux` and
	`fastFluence` are not created.
	
	\note
	The abscence of the `fastFlux` and `fastFluence` fields might create conflicts
	with other models that partially or entirely depends on them (e.g. some of the
	material properties models). If you want to use those models bust still neglect 
    the fast flux and fluence evolution over time, use the `constant` fastFlux model.
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for selecting the `none`
    fast flux model
	
	\verbatim
	    fastFlux none;
	\endverbatim

SourceFiles
    fastFlux.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021


\*---------------------------------------------------------------------------*/

#ifndef flux_H
#define flux_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "materials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class fastFlux Declaration
\*---------------------------------------------------------------------------*/

class fastFlux
{
    // Private data
    
        //-Private Member Functions

        //- Disallow default bitwise copy construct
        fastFlux(const fastFlux&);

        //- Disallow default bitwise assignment
        void operator=(const fastFlux&);
        
protected:  
    
    // Protected Data
        
        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to the materials class
        const materials& mat_;

        //- FastFluxOptions dictionary
        const dictionary fastFluxDict_;

    // Protected member functions
        
        //- Update fluence
        virtual void advanceFluence()
        {
            NotImplemented;
        };

public:

    //- Runtime type information
    TypeName("none");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fastFlux,
            dictionary,
            ( 
                const fvMesh& mesh,
                const materials& materials,
                const dictionary& fastFluxDict
            ),
            (mesh, materials, fastFluxDict)
        );

    // Constructors

        //- Construct from mesh, materials and dict
        fastFlux
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& fastFluxDict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<fastFlux> New
        (        
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& solverDict
        );


    //- Destructor
    virtual ~fastFlux();


    // Access functions
      
        //- Return mesh reference (not implemented)
        const fvMesh& mesh() const
        {
            return mesh_;
        } 
    
        //- Return fast flux field reference (not implemented)
        const volScalarField& flux() const
        {
            NotImplemented;
            return volScalarField::null();
        };
    
        //- Return fast fluence field reference
        const volScalarField& fastFluence() const
        {
            NotImplemented;
            return volScalarField::null();
        };

    // Member Functions

        //- Do nothing
        virtual void correct()
        {};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
