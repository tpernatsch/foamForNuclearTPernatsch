/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::virtualMass

Description
    Class to handle fluid-fluid interfacial area

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    virtualMass.C
    newInterfacialAreaModel.C

\*---------------------------------------------------------------------------*/

#ifndef virtualMass_H
#define virtualMass_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "fvMatrix.H"
#include "virtualMassCoefficientModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class FFPair;

/*---------------------------------------------------------------------------*\
            Class virtualMass Declaration
\*---------------------------------------------------------------------------*/

class virtualMass
:
    public IOdictionary
{
protected:

    const fvMesh& mesh_;

    const FFPair& pair_;

    const fluid& fluid1_;

    const fluid& fluid2_;

    const volVectorField& U1_;
    
    const volVectorField& U2_;
    
    const surfaceScalarField& phi1_;
    
    const surfaceScalarField& phi2_;

    autoPtr<virtualMassCoefficientModel> VmPtr_;

    volScalarField Vm_;

    typedef HashTable<fvVectorMatrix, word, word::hash> VmForceTable;

    VmForceTable VmForces_;

public:

    //- Runtime type information
    TypeName("virtualMass");

    // Constructors

        virtualMass
        (
            const FFPair& pair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~virtualMass(){}

    // Member Functions

        //-
        virtual void correct();

        //- Access

            const fvVectorMatrix& force(const word& UName) const
            {
                return VmForces_[UName];
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
