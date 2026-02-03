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
    Foam::HastelloyN

Description
    HastelloyN thermomechanical properties are taken from different sources.

    In general, a certain material property or behavioral model can 
    be selected with the specific keyword. Otherwise, default models will be used.

    Nevertheless, for HN the default "constant" model is used for density,
    heat capacity, emissivity and Poisson Ratio. Thus the user needs
    to introduce their values. Typical values are shown in the example below.

Usage
    In solverDict file:
    \verbatim
    materials
    {
        //- Name of cellZone with HastelloyN material model
        cladding
        {
            material                HastelloyN;
            Tref                    Tref [ 0 0 0 1 0 ] 293;

            // For the following properties the "constant" model is 
            // selected by default, but it can be changed if needed
            rho         rho     [1 -3 0 0 0]    8860;
            Cp          Cp      [0 2 -2 -1 0]   420;
            emissivity emissivity [0 0 0 0 0]   0.25;        
            nu          nu      [0 0 0 0 0]     0.3;  

            // The following models are already selected by default
            // but can be changed if needed
            conductivityModel SwindemanHastelloyN;
            YoungModulusModel WatrousHastelloyN;
            thermalExpansionModel SwindemanHastelloyN;
            swellingModel WrightShamHastelloyN;
            
            rheologyModel               elasticity;   
            
        }
 
        ...

    }
    \endverbatim

SourceFiles
    HastelloyN.C

 \mainauthor
    Alejandra de Lara - University of Cambridge

 \contribution
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland,
    Laboratory for Reactor Physics and Systems Behaviour)\n

 \date
     Sep 2023

\*---------------------------------------------------------------------------*/

#ifndef HastelloyN_H
#define HastelloyN_H

#include "materialModel.H"
#include "volFields.H"

#include "swellingModel.H"
#include "phaseTransitionModel.H"

#include "failureModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class HastelloyN Declaration
\*---------------------------------------------------------------------------*/

class HastelloyN
:
    public materialModel
{
    // Private data

        //- irradiationGrowth model
        autoPtr<swellingModel> swelling_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        HastelloyN(const HastelloyN&);

        //- Disallow default bitwise assignment
        void operator=(const HastelloyN&);

        
    //- Scaling parameters for UQ and SA
    
        //- Group name for parameter definition
        declareNameWithVirtual(HastelloyN, group_, group);

        //- Density scaling factor [-]
        declareParameterWithVirtual(HastelloyN, F_rho_, F_rho);
        
        //- Specific heat scaling factor [-]
        declareParameterWithVirtual(HastelloyN, F_Cp_, F_Cp);
    
        //- Thermal conductivity scaling factor [-]
        declareParameterWithVirtual(HastelloyN, F_k_, F_k);
    
        //- Emissivity scaling factor [-]
        declareParameterWithVirtual(HastelloyN, F_emissivity_, F_emissivity);
    
        //- Elastic Young's modulus scaling factor [-]
        declareParameterWithVirtual(HastelloyN, F_E_, F_E);
    
        //- Poisson's ratio scaling factor [-]
        declareParameterWithVirtual(HastelloyN, F_nu_, F_nu);
        
        //- Themral expansion (alpha*DeltaT) scaling factor  [-]
        declareParameterWithVirtual(HastelloyN, F_alphaT_, F_alphaT);

        // //- Density offset [kg/m3]
        declareParameterWithVirtual(HastelloyN, delta_rho_, delta_rho);
        
        //- Specific heat offset [J/kg.K]
        declareParameterWithVirtual(HastelloyN, delta_Cp_, delta_Cp);
    
        //- Thermal conductivity offset [W/m.K]
        declareParameterWithVirtual(HastelloyN, delta_k_, delta_k);
    
        //- Emissivity offset [-]
        declareParameterWithVirtual(HastelloyN, delta_emissivity_, delta_emissivity);
    
        //- Elastic Young's modulus offset [N/m2]
        declareParameterWithVirtual(HastelloyN, delta_E_, delta_E);
    
        //- Poisson's ratio offset [-]
        declareParameterWithVirtual(HastelloyN, delta_nu_, delta_nu);
        
        //- Themral expansion (alpha*DeltaT) offset [-]
        declareParameterWithVirtual(HastelloyN, delta_alphaT_, delta_alphaT);


public:

    //- Runtime type information
    TypeName("HastelloyN");


    // Constructors

        //- Construct from mesh, global options, materials, and dictionary
        HastelloyN
        (
            const fvMesh& mesh,
            const dictionary& materialModelDict,
            const labelList& addr
        );


    //- Destructor
    ~HastelloyN();


    // Access Functions


    // Member Functions
    
        //- Update the behavioral models
        virtual void correctBehavioralModels
        (
            const scalarField& T
        );
        
        //- Return the tensor of additional strains
        virtual void additionalStrains( symmTensorField& sf ) const;

        //- Return name of the material
        virtual word name() const
        {
            return "HastelloyN";
        };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

