/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::burnup

Description
	As the name implies, the `none` burnup class in OFFBEAT allows you to switch off
	the burnup model and have a simulation where burnup is not at all considered.
	This means that the burnup field `Bu` is not created nor present in the simulation.
	
	\warning
	The abscence of the `Bu` field might create conflicts with other models that
	partially or entirely depends on the burnup (e.g. some fuel conductivity models).
	If you want to use those models bust still neglect the burnup evolution over time, 
    use the `constant` burnup.
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for activating the
	`none` burnup model
	
	\verbatim
	    burnup none;
	\endverbatim
	
SourceFiles
    burnup.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef burnup_H
#define burnup_H

#include "primitiveFields.H"
#include "volFields.H"
#include "materials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class burnup Declaration
\*---------------------------------------------------------------------------*/

class burnup
{
    // Private data
    
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        burnup(const burnup&);

        //- Disallow default bitwise assignment
        void operator=(const burnup&);

protected:
        
    //- Reference to mesh
    const fvMesh& mesh_;

    //- Reference to the materials class
    const materials& mat_;

    //- Input/options dictionary
    const dictionary burnupDict_;

public:

    //- Runtime type information
    TypeName("none");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        burnup,
        dictionary,
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& burnupDict
        ),
        (mesh, mat, burnupDict)
    );

    // Constructors

        //- Construct from mesh, materials and dictionary
        burnup
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& burnupDict
        );

    // Selectors

        //- Select from dictionary
        static autoPtr<burnup> New
        (
            const fvMesh& mesh,
            const materials& mat,
            const dictionary& solverDict
        );


        //- Destructor
        virtual ~burnup();


    // Member Functions

    
        //- Update the burnup   
        virtual void correct()
        {};

        //- Return a constant reference to the burnup field
        virtual const volScalarField& Bu() const
        {
            NotImplemented

            return volScalarField::null();
        };

        //- Return the SigmaA_ 
        virtual const scalarField SigmaA() const
        {
            NotImplemented

            return volScalarField::null();
        };

        //- Return the diffusion coeff field 
        virtual const scalarField neutronD() const
        {
            NotImplemented

            return volScalarField::null();
        };
        
        //- Return next time increment according to burnup model
        virtual scalar nextDeltaT()
        {
            return GREAT;
        };

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
