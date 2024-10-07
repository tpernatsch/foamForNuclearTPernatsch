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
    Foam::inconel600

Description
    Inconel600 thermoMechanical properties, taken mostly from the web.

SourceFiles
    inconel600.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef inconel600_H
#define inconel600_H

#include "materialModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class inconel600 Declaration
\*---------------------------------------------------------------------------*/

class inconel600
:
    public materialModel
{
    // Private data
    
        //- Constant density [kg/m3]
        dimensionedScalar rho_;
    
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        inconel600(const inconel600&);

        //- Disallow default bitwise assignment
        void operator=(const inconel600&);

public:

    //- Runtime type information
    TypeName("inconel600");


    // Constructors
    
        //- Construct from mesh, materials, and dictionary
        inconel600
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~inconel600();


    // Access Functions


    // Member Functions

        //- Modifies the density [kg/m3] as a function of temperature for the supplied
        //  addressing
        virtual void rho
        (scalarField& sf, const scalarField& T, const labelList& addr);
        
        //- Modifies the specific heat capacity [J/kg.K] as a function of temperature
        //  for the supplied addressing
        virtual void Cp
        (scalarField& sf, const scalarField& T, const labelList& addr);
        
        //- Modifies the thermal conductivity [W/m.K] as a function of temperature for 
        //  the supplied addressing
        virtual void k
        (scalarField& sf, const scalarField& T, const labelList& addr);
        
        //- Modifies the emissivity [-] as a function of temperature for the supplied
        //  addressing
        virtual void emissivity
        (scalarField& sf, const scalarField& T, const labelList& addr);
        
        //- Modifies the elastic modulus [Pa] as a function of temperature for the
        //  supplied addressing
        virtual void E
        (scalarField& sf, const scalarField& T, const labelList& addr);
        
        //- Modifies the Poisson's ratio [-] as a function of temperature for the
        //  supplied addressing
        virtual void nu
        (scalarField& sf, const scalarField& T, const labelList& addr);
        
        //- Modifies the linear expansion dL/L [-] as a function of temperature
        //  for the supplied addressing
        virtual void alphaT
        (symmTensorField& sf, const scalarField& T, const labelList& addr) const;

        //- True if material is of "fuel" type
        virtual word name() const
        {
            return "inconel600";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
