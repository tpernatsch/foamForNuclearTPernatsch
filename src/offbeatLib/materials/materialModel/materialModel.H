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
    Foam::materialModel

Description
    Parent class for all materials. Some of its functions are pure virtual and
    must be implemented in the daughter classes.

SourceFiles
    materialModel.C

\todo
    - find better way to read isotopic compositions --> maybe create mother 
    class "fuelMaterialModel.C" and from there create Uo2 and PuUO2

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

    

\*---------------------------------------------------------------------------*/

#ifndef materialModel_H
#define materialModel_H

#include "primitiveFields.H"
#include "objectRegistry.H"
#include "volFields.H"

#include "densityModel.H"
#include "heatCapacityModel.H"
#include "conductivityModel.H"
#include "emissivityModel.H"
#include "YoungModulusModel.H"
#include "PoissonRatioModel.H"
#include "thermalExpansionModel.H"
#include "damageModel.H"

#include "userParameters.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class materialModel Declaration
\*---------------------------------------------------------------------------*/

class materialModel
{
    // Private data
    
    // Private Member Functions

    //- Disallow default bitwise copy construct
    materialModel(const materialModel&);

    //- Disallow default bitwise assignment
    void operator=(const materialModel&);

protected:
    
    // Protected data
    
    //- Reference to mesh
    const fvMesh& mesh_;

    //- Material model dictionary
    const dictionary materialModelDict_;

    //- Ref to address of the material model
    const labelList& addr_;

    //- density model
    autoPtr<densityModel> density_;

    //- heat Capacity model
    autoPtr<heatCapacityModel> heatCapacity_;
    
    //- conductivity model
    autoPtr<conductivityModel> conductivity_;
    
    //- emissivity model
    autoPtr<emissivityModel> emissivity_;
    
    //- YoungModulus model
    autoPtr<YoungModulusModel> YoungModulus_;
    
    //- PoissonRatio model
    autoPtr<PoissonRatioModel> PoissonRatio_;
    
    //- thermal Expansion model
    autoPtr<thermalExpansionModel> thermalExpansion_;
    
    //- damage or softening model
    autoPtr<damageModel> damageModel_;

    //- Scaling parameters for UQ and SA
        
        //- Group name for parameter definition
        declareNameWithVirtual(materialModel, group_, group);
    
        //- Density scaling factor [-]
        declareParameterWithVirtual(materialModel, F_rho_, F_rho);
        
        //- Specific heat scaling factor [-]
        declareParameterWithVirtual(materialModel, F_Cp_, F_Cp);
    
        //- Thermal conductivity scaling factor [-]
        declareParameterWithVirtual(materialModel, F_k_, F_k);
    
        //- Emissivity scaling factor [-]
        declareParameterWithVirtual(materialModel, F_emissivity_, F_emissivity);
    
        //- Elastic Young's modulus scaling factor [-]
        declareParameterWithVirtual(materialModel, F_E_, F_E);
    
        //- Poisson's ratio scaling factor [-]
        declareParameterWithVirtual(materialModel, F_nu_, F_nu);
        
        //- Themral expansion (alpha*DeltaT) scaling factor  [-]
        declareParameterWithVirtual(materialModel, F_alphaT_, F_alphaT);

        // //- Density offset [kg/m3]
        declareParameterWithVirtual(materialModel, delta_rho_, delta_rho);
        
        //- Specific heat offset [J/kg.K]
        declareParameterWithVirtual(materialModel, delta_Cp_, delta_Cp);
    
        //- Thermal conductivity offset [W/m.K]
        declareParameterWithVirtual(materialModel, delta_k_, delta_k);
    
        //- Emissivity offset [-]
        declareParameterWithVirtual(materialModel, delta_emissivity_, delta_emissivity);
    
        //- Elastic Young's modulus offset [N/m2]
        declareParameterWithVirtual(materialModel, delta_E_, delta_E);
    
        //- Poisson's ratio offset [-]
        declareParameterWithVirtual(materialModel, delta_nu_, delta_nu);
        
        //- Themral expansion (alpha*DeltaT) offset [-]
        declareParameterWithVirtual(materialModel, delta_alphaT_, delta_alphaT);

public:

    //- Runtime type information
    TypeName("materialModel");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        materialModel,
        dictionary,
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        ),
        (mesh, materialModelDict, addr)
    );

    // Constructors

    //- Construct from mesh, materials, and address
    materialModel
    (
        const fvMesh& mesh, 
        const dictionary& materialModelDict,
        const labelList& addr
    );


    // Selectors

    //- Select from dictionary
    static autoPtr<materialModel> New
    (
        const fvMesh& mesh, 
        const dictionary& materialModelDict,
        const labelList& addr
    );

    //- Destructor
    virtual ~materialModel();

    // Access Functions
    virtual const dictionary materialModelDict() const
    {
        return materialModelDict_;
    }

    // Member Functions - ThermoMech Properties

    //- Return the density [kg/m3] as a function of temperature for the supplied
    //- addressing
    virtual void rho
    (scalarField& sf, const scalarField& T);
    
    //- Return the specific heat capacity [J/kg.K] as a function of temperature
    //- for the supplied addressing
    virtual void Cp
    (scalarField& sf, const scalarField& T);
    
    //- Return the thermal conductivity [W/m.K] as a function of temperature for 
    //- the supplied addressing
    virtual void k
    (scalarField& sf, const scalarField& T);
    
    //- Return the emissivity [-] as a function of temperature for the supplied
    //- addressing
    virtual void emissivity
    (scalarField& sf, const scalarField& T);
    
    //- Return the elastic modulus [Pa] as a function of temperature for the
    //- supplied addressing
    virtual void E
    (scalarField& sf, const scalarField& T);
    
    //- Return the Poisson's ratio [-] as a function of temperature for the
    //- supplied addressing
    virtual void nu
    (scalarField& sf, const scalarField& T);
    
    //- Return the linear expansion dL/L [-] as a function of temperature
    //- for the supplied addressing
    virtual void alphaT
    (symmTensorField& sf, const scalarField& T, const labelList& addr) const;
    
    //- Softens stiffness E and Poisson Ratio nu depending on damage model
    virtual void soften
    (scalarField& E, scalarField& nu);
    

    virtual void correctDamage
    (
        volScalarField& sigmaHyd, 
        volSymmTensorField& sigmaDev,
        volSymmTensorField& epsilonEl
    );

    //- Apply user parameters to given field
    template<class T>
    void applyParameters
    (
        Field<T>& sf, 
        const labelList& addr, 
        const dimensionedScalar& F, 
        const dimensionedScalar& delta
    ) const
    {    
        forAll(addr, i)
        {
            label cellI(addr[i]);

            sf[cellI] = F.value()*sf[cellI] + delta.value()*pTraits<T>::one;
        }
    };

    // Member Functions - Behavioral Models
    
    //- Update the tensor of additional strains
    virtual void correctBehavioralModels( const scalarField& T ){};

    //- Check if failure occurred in the material
    virtual void checkFailure(){};
    
    //- Return the tensor of additional strains
    virtual void additionalStrains( symmTensorField& sf ) const{};

    //- Return material name
    virtual word name() const {return "materialModel";};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
