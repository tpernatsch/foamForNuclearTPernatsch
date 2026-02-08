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
    Foam::UPuO2

Description
    UPuO2 material class. 

    The thermomechanical properties are taken from different sources, as Bison 
    manual, Frapcon manual or Relap manual.

    The material property and behavioral models can selected with the specific
    keywords. Otherwise, default models will be used.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        //- Name of cellZone with UPuO2 material model
        fuel
        {
            material UPuO2;

            densityModel            UPuO2Constant;
            heatCapacityModel       UPuO2Matpro;
            conductivityModel       UPuO2;
            emissivityModel         UO2Relap;
            YoungModulusModel       UPuO2Matpro;
            PoissonRatioModel       UPuO2Constant;
            thermalExpansionModel   UPuO2Matpro;

            densificationModel      UO2Frapcon;
            swellingModel           UO2Frapcon;
            relocationModel         UO2Frapcon;
            failureModel            none;
        }
        ...

    }
    \endverbatim

SourceFiles
    UPuO2.C

\todo

    - No correlation for (U,Pu)O2 emissivity is available. It is therefore 
    recommended to use the one for UO2 fuel. Source:
    (https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-19417.pdf).
    - Based on the limited available in-reactor data for MOX, the same 
    correlations for UO2 and UPuO2 are used for densification and swelling.
    Source: Matprov11
    - The same goes for the relocation, no model for MOX is available
    and therefore the UO2 model is used.

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef UPuO2_H
#define UPuO2_H

#include "fuelMaterial.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class UPuO2 Declaration
\*---------------------------------------------------------------------------*/

class UPuO2
:
    public fuelMaterial
{
    // Private data

    // Private Member Functions

        //- Initial Pu/HM ratio
        const scalar PuMetalRatio_;

        //- Initial U/HM ratio
        const scalar UMetalRatio_;

        //- Disallow default bitwise copy construct
        UPuO2(const UPuO2&);

        //- Disallow default bitwise assignment
        void operator=(const UPuO2&);

        
    //- Scaling parameters for UQ and SA
    
        //- Group name for parameter definition
        declareNameWithVirtual(UPuO2, group_, group);

        //- Density scaling factor [-]
        declareParameterWithVirtual(UPuO2, F_rho_, F_rho);
        
        //- Specific heat scaling factor [-]
        declareParameterWithVirtual(UPuO2, F_Cp_, F_Cp);
    
        //- Thermal conductivity scaling factor [-]
        declareParameterWithVirtual(UPuO2, F_k_, F_k);
    
        //- Emissivity scaling factor [-]
        declareParameterWithVirtual(UPuO2, F_emissivity_, F_emissivity);
    
        //- Elastic Young's modulus scaling factor [-]
        declareParameterWithVirtual(UPuO2, F_E_, F_E);
    
        //- Poisson's ratio scaling factor [-]
        declareParameterWithVirtual(UPuO2, F_nu_, F_nu);
        
        //- Themral expansion (alpha*DeltaT) scaling factor  [-]
        declareParameterWithVirtual(UPuO2, F_alphaT_, F_alphaT);

        // //- Density offset [kg/m3]
        declareParameterWithVirtual(UPuO2, delta_rho_, delta_rho);
        
        //- Specific heat offset [J/kg.K]
        declareParameterWithVirtual(UPuO2, delta_Cp_, delta_Cp);
    
        //- Thermal conductivity offset [W/m.K]
        declareParameterWithVirtual(UPuO2, delta_k_, delta_k);
    
        //- Emissivity offset [-]
        declareParameterWithVirtual(UPuO2, delta_emissivity_, delta_emissivity);
    
        //- Elastic Young's modulus offset [N/m2]
        declareParameterWithVirtual(UPuO2, delta_E_, delta_E);
    
        //- Poisson's ratio offset [-]
        declareParameterWithVirtual(UPuO2, delta_nu_, delta_nu);
        
        //- Thermal expansion (alpha*DeltaT) offset [-]
        declareParameterWithVirtual(UPuO2, delta_alphaT_, delta_alphaT);

public:

    //- Runtime type information
    TypeName("UPuO2");


    // Constructors

        //- Construct from mesh, materials, and dictionary
        UPuO2
        (
            const fvMesh& mesh, 
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~UPuO2();


    // Access Functions

        //- Return initial Pu/HM ratio
        virtual const scalar PuMetalRatio() const
        {
            return PuMetalRatio_;
        };

        //- Return initial U/HM ratio
        virtual const scalar UMetalRatio() const
        {
            return UMetalRatio_;
        };

    // Member Functions
    
        //- Update the behavioral models
        virtual void correctBehavioralModels( const scalarField& T );

        //- Check if failure occurred 
        virtual void checkFailure();
        
        //- Return the tensor of additional strains
        virtual void additionalStrains( symmTensorField& sf ) const;

        //- Return name of the material
        virtual word name() const
        {
            return "UPuO2";
        }; 
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
