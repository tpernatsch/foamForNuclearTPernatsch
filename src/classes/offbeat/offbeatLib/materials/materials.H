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
    Foam::materials

Description
    Parent class that handles the evolution of material properties and 
    behavioral models. 
    
    The material models are stored in a PtrList. The addressing of each material
    is stored in a labelListList.
    
    The mother class is not runtime-selectable.

SourceFiles
    materials.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef materials_H
#define materials_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOdictionary.H"

#include "primitiveFields.H"
#include "materialModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class materials Declaration
\*---------------------------------------------------------------------------*/
class materialModel;

class materials
{
    // Private data

    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        materials(const materials&);

        //- Disallow default bitwise assignment
        void operator=(const materials&);

protected:
    
    // Protected data   

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Material dictionary
        const dictionary materialsDict_;

        //- List of materials
        PtrList<materialModel> matList_;

        //- List of material addressing
        mutable labelListList matAddrList_;
        
        //- Temperature pointer
        const volScalarField* T_;
        
        //- Density [kg/m3]
        volScalarField rho_;
        
        //- Specific heat capacity [J/kg.K]
        volScalarField Cp_;
        
        //- Thermal conductivity [W/m.K]
        volScalarField k_;
        
        //- emissivity [-]
        volScalarField emissivity_;
        
        //- Elastic (Young's) modulus [Pa]
        volScalarField E_;
        
        //- Poisson's ratio [-]
        volScalarField nu_;

        //- First Lamé parameter modulus [Pa]
        volScalarField mu_;
        
        //- Second Lamé parameter modulus [Pa]
        volScalarField lambda_;
        
        //- Bulk modulus times 3 [Pa]
        volScalarField threeK_;

        //- Switch to write all the material thermomechanical properties
        Switch writeProperties_;
        
        //- Maximum material properties increase for adjustable time step
        scalar maxDeltaMatTol_;

    // Protected member functions

        //- Update addressing for each material
        virtual void calcAddressing(const dictionary& materialsDict) const = 0;

public:

    //- Runtime type information
    TypeName("materials");
    

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            materials,
            dictionary,
            (
                const fvMesh& mesh, 
                const dictionary& materialsDict
            ),
            (mesh, materialsDict)
        );


    // Constructors

        //- Construct from material dictionary
        materials
        (
            const fvMesh& mesh,
            const dictionary& materialsDict
        );

        //- Return a pointer to a new materials
        static autoPtr<materials> New
        (
            const fvMesh& mesh,
            const dictionary& solverDict
        );


    //- Destructor
    virtual ~materials();


    // Access Functions

        //- Return a constant reference to the density [kg/m3]
        inline const volScalarField& rho() const
        {
            return rho_;
        }
        
        //- Return a constant reference to the specific heat [J/kg.K]
        inline const volScalarField& Cp() const
        {
            return Cp_;
        }
        
        //- Return a constant reference to the thermal conductivity [W/m.K]
        inline const volScalarField& k() const
        {
            return k_;
        }

        //- Return a constant reference to the emissivity [-]
        inline const volScalarField& emissivity() const
        {
            return emissivity_;
        }

        //- Return a constant reference to the elastic modulus [Pa]
        inline const volScalarField& E() const
        {
            return E_;
        }
        
        //- Return a constant reference to the Poisson's ratio [-]
        inline const volScalarField& nu() const
        {
            return nu_;
        }
        
        //- Lame's first parameter
        inline const volScalarField& lambda() const
        {
            return lambda_;
        }
        
        //- Lame's second parameter (shear modulus)
        inline const volScalarField& mu() const
        {
            return mu_;
        }
        
        //- Bulk modulus
        inline const volScalarField& threeK() const
        {
            return threeK_;
        }
        
        //- Return reference to materials list
        inline const PtrList<materialModel>& materialsList() const
        {
            return matList_;
        }
        
        //- Return reference to material addressing list
        virtual const labelListList& matAddrList() const 
        {
            // Topology has changed
            if(mesh_.topoChanging())
            {
                calcAddressing(materialsDict_);
            }

            return matAddrList_;
        };    


    // Member Functions

        //- Correct/update the thermoMechanical properties
        virtual void correctThermoMechProperties() = 0;

        //- Correct/update the behavioral models
        virtual void correctBehavioralModels() = 0;

        //- Check material failure occurrence
        virtual void checkFailure() = 0;
        
        //- Linear thermal expansion dL/L
        virtual tmp<volSymmTensorField> thermalExpansion() const = 0;
        
        //- Additional strain components per material
        virtual tmp<volSymmTensorField> additionalStrains() const = 0; 

        //- Update density of the cells the mesh is updated
        virtual void updateDensity();        

        //- Calculate max deltaT based on material properties relative change
        virtual scalar nextDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
