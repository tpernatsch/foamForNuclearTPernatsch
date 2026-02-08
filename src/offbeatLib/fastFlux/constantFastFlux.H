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
    Foam::constantFastFlux

Description
	The `constant` fastFlux class in OFFBEAT allows you to set the fast flux
	and fast fluence fields in files located the starting time folder. If the
	`fastFlux` and `fastFluence` files are present in the starting time folder, the
	internal field and boundary conditions are set accordingly. Otherwise, the
	fields are set to 0 <code>n/cm<sup>2</sup>/s</code>  and 0 <code>n/cm<sup>2</sup></code> respectively
	in each cell, with zeroGradient BCs in all non-empty/-wedge patches. Note that
	when selecting this fast flux model, the fast flux field will remain the same as
	it was defined in the starting time folder (or set by default).
	
	\note
	Even when `constant` is selected, the fast fluence keeps being updated at
	each time step as:
	\f$fastFluence = fastFluence.oldTime() + fastFlux \cdot \Delta t\f$

    \warning
    For backward-compatibility, `constant` can be still used as the typename
    for the `constant` fast flux model. However, its use is deprecated.
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for selecting the `constant`
    fast flux model
	
	\verbatim
	    fastFlux constant;
	\endverbatim

SourceFiles
    constantFastFlux.C

\todo
    - material zones provided or decided based on material.isFuel() or everywhere
    - do we need also thermalFlux classes?
    - do we need the profile to be normalized?
    - do we need a radial profile for the fast/thermal flux?    

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef constantFastFlux_H
#define constantFastFlux_H

#include "fastFlux.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantFastFlux Declaration
\*---------------------------------------------------------------------------*/

class constantFastFlux
:
    public fastFlux
{ 
    // Private data
    
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        constantFastFlux(const constantFastFlux&);

        //- Disallow default bitwise assignment
        void operator=(const constantFastFlux&);

protected:
    
    // Protected Data    
        
        //- FastFluence in n/cm2
        volScalarField fastFluence_;
        
        //- FastFlux in n/cm2/s
        volScalarField flux_;    
            
        //- Value of current time
        scalar currentTime_;   

    // Protected member functions
        
        //- Update fluence
        virtual void advanceFluence();

public:

    //- Runtime type information
    TypeName("constant");

    // Constructors

        //- Construct from mesh, materials and dict
        constantFastFlux
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& fastFluxDict
        );


    //- Destructor
    ~constantFastFlux();


    // Access functions
    
        //- Return fast flux field reference
        const volScalarField& flux() const
        {
           return flux_;
        };
    
        //- Return fast fluence field reference
        const volScalarField& fastFluence() const
        {
           return fastFluence_;
        };
        
    // Member Functions

        //- Update fast flux and fluence
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
