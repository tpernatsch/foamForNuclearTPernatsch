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
    Foam::materialsByZone

Description
    The materials are identified with the cellZones. However, the same material
    model (e.g. UO2, MOX, Zircaloy etc) can be used in different cellZones.

    The material properties are provided in a subDictionary called "materials".
    The class will check that all cellZones are assigned with one (and only one)
    material model.

SourceFiles
    materialsByZone.C

Usage
    In solverDict file:
    \verbatim
    materialProperties byZone;

    // List of materials, one per cellZone.
    materials
    {
        fuel
        {
            material UO2;
            ...
        }

        cladding
        {
            material Zircaloy;
            ...
        }
    }
    \endverbatim

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef materialsByZone_H
#define materialsByZone_H

#include "materials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class materialsByZone Declaration
\*---------------------------------------------------------------------------*/

class materialsByZone
:
    public materials
{
    // Private data
    
    
    // Private Member Functions
        
        //- Check that all cells are accounted for in zones
        void checkZoneDefinition();

        //- Disallow default bitwise copy construct
        materialsByZone
        (
            const materialsByZone&
        );

        //- Disallow default bitwise assignment
        void operator=(const materialsByZone&);

protected:

    // Protected member functions

        //- Update addressing for each material
        virtual void calcAddressing(const dictionary& materialsDict) const;

public:

    //- Runtime type information
    TypeName("byZone");

    
    // Constructors

        //- Construct from mesh and dictionary
        materialsByZone
        (
            const fvMesh& mesh,
            const dictionary& materialsDict
        );

    
    //- Destructor
    ~materialsByZone();


    // Member Functions

        //- Correct/update the properties
        virtual void correctThermoMechProperties();

        //- Correct/update the behavioral models
        virtual void correctBehavioralModels();
        
        //- Check material failure occurrence
        virtual void checkFailure();
        
        //- Linear thermal expansion dL/L
        virtual tmp<volSymmTensorField> thermalExpansion() const; 
        
        //- Additional strain components per material
        virtual tmp<volSymmTensorField> additionalStrains() const; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
