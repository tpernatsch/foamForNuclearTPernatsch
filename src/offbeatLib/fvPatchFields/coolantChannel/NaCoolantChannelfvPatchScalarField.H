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
    Foam::NaCoolantChannelfvPatchScalarField

Description
    Coolant subchannel module to provide thermal boundary condition for Na 
    cooled fuel rods. The model solves the 1D energy equation in the coolant 
    along the axial direction to calculate the coolant bulk temperature axial
    profile. The surface temperature on the cladding is computed through the 
    heat transfer coefficient \f$h\f$. To determine \f$h\f$, different Nusselt
    correlations are implemented in the model.

________________________________________________________________________________

\par Options
    The `NaCoolantChannel` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    Parameters in the patch subdictionary for `NaCoolantChannel`:

    - <b>`flowArea`</b> - Flow area in m\f$^2$.
    - <b>`hydraulicDiameter`</b> - Hydraulic diameter in m.
    - <b>`massFlowRate`</b> - Mass flowrate in kg/s (time-dependent list is also possible).
    - <b>`inletTemperature`</b> - Coolant inlet temperature (K) (time-dependent list is also possible).
    - <b>`NusseltModel`</b> - Selected model for Nu computation.
    - <b>`conductivityModel`</b> - Selected model for \f$k_{Na}$ computation.
    - <b>`viscosityModel`</b> - Selected model for \f$\nu$ computation.
    - <b>`heatCapacityModel`</b> - Selected model for \f$\C_p$ computation.

________________________________________________________________________________

Usage
    In T boundary condition dictionary:
    \verbatim
        cladOuterSurface
        {
            type                    NaCoolantChannel;

            // Flow area in m2
            flowArea                20e-6;
            
            // Hydraulic diameter in m
            hydraulicDiameter       4.0e-3;
            
            // Mass flowrate in kg/s (time-dependent list is also possible)
            massFlowRate            0.15;
            
            // Coolant inlet temperature (K) (time-dependent list is also possible)
            inletTemperature        565;

            // Under-relaxation factor for heat transfer coefficient
            relaxHTC                1.0;

            // Nusselt correlations (Default: Schad)
            NusseltModel            Schad;
            // NusseltModel            Lyon;
            // NusseltModel            SebanShimazaki;
            // NusseltModel            DittusBoelter;
            // NusseltModel            constant;
            
            // Only if NusseltModel Schad is selected
            pitchOverDiameter       1.2;
            
            // Sodium Conductivity Model (Default : FinkLeibowitz) 
            conductivityModel       FinkLeibowitz;
            
            // Sodium Heat Capacity Model (Default : Vargaftic) 
            heatCapacityModel       Vargaftic;

            // Sodium viscosity model (default : Shpilrain)
            viscosityModel          Shpilrain; 

            // Switch to write Tcoolant and enthalpy volFields
            writeCoolantT           true; 

            // Initialize cladding temperature
            value                   uniform 565;
        }
    \endverbatim 

SourceFiles
    NaCoolantChannelfvPatchScalarField.C

\mainauthor
    E. Brunetto, A. Sz. Vanyi - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    July 2024

\*---------------------------------------------------------------------------*/

#ifndef NaCoolantChannelfvPatchScalarField_H
#define NaCoolantChannelfvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "axialProfile.H"
#include "InterpolateTables.H"
#ifdef OPENFOAMESI
    #include "Function1.H"
#endif     
#include "Table.H"
#include "sliceMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class NaCoolantChannelfvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class NaCoolantChannelfvPatchScalarField
:
    public fvPatchField<scalar>
{
    // Private data
        
        //- Name of thermal conductivity field
        word kappaName_;
        
        //- T0 patch profile @ current timestep - working variable
        scalarField T0_;

        //- h patch profile @ current timestep - for write()
        scalarField h_;
        
        //- Boundary heat transfer factor - working variable
        scalarField alpha_;

        //- Subchannel flow area
        scalar pOverD_;

        //- Subchannel flow area
        scalar A_;

        //- Relaxation factor of the HTC
        scalar relaxHTC_;

        //- Hydraulic diameter
        scalar DH_;

        //- Coolant mass flow rate
        scalar massFlowRate_;

#ifdef OPENFOAMFOUNDATION
        //- Coolant mass flow rate time-dependent list
        autoPtr<Function1s::Table<scalar>> massFlowRateList_;
#elif OPENFOAMESI
        //- Coolant mass flow rate time-dependent list
        autoPtr<Function1<scalar>> massFlowRateList_;
#endif     

        //- Inlet temperature
        scalar inletTemperature_;

#ifdef OPENFOAMFOUNDATION
        //- Inlet temperature time-dependent list
        autoPtr<Function1s::Table<scalar>> inletTemperatureList_;
#elif OPENFOAMESI
        //- Inlet temperature time-dependent list
        autoPtr<Function1<scalar>> inletTemperatureList_;
#endif 

        //- Correlation for liquid sodium conductivity
        word modelKappa_;

        //- Constant Kappa value
        scalar constKappa_;
        
        //- Correlation for liquid sodium specific heat capacity
        word modelCP_;

        //- Constant cp value
        scalar constCp_;
        
        //- Correlation for liquid sodium viscosity
        word modelViscosity_;

        //- Constant viscosity value
        scalar constVisc_;

        //- Correlation for Nusselt number
        word modelNusselt_;

        //- Constant Nusselt value
        scalar constNusselt_;

        //- Switch to write T coolant as a field (for postprocessing)
        Switch writeCoolantT_;

        //- Autoptr to T coolant field (for postprocessing)
        autoPtr<volScalarField> Tcoolant_;
        
        //- Autoptr to specific enthalpy coolant field (for postprocessing)
        autoPtr<volScalarField> Hcoolant_;

        //- Heat flux scalar field (for debug)
        scalarField heatFlux_;

        //- Reference to slice mapper
        const sliceMapper* mapper_;   

    // Private functions

        // Copmute slice heights
        virtual scalarField computeSliceHeights();

        // Returns conductivity of liquid Sodium at a given Temperature
        virtual scalar NaConductivity(scalar T);

        // Returns enhalpy of liquid Sodium at a given Temperature
        virtual scalar NaEnthalpy(scalar T);

        // Returns Specific heat capacity of liquid Sodium at a given Temperature
        virtual scalar NaCp(scalar T);

        // Returns dynamic viscosity of liquid Sodium at a given Temperature
        virtual scalar NaVisc(scalar T);

        // Returns Pr number of liquid Sodium at a given Temperature
        virtual scalar NaPrandtl(scalar T);
        
public:

    //- Runtime type information
    TypeName("NaCoolantChannel");


    // Constructors

        //- Construct from patch and internal field
        NaCoolantChannelfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        NaCoolantChannelfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given NaCoolantChannelfvPatchScalarField
        //  onto a new patch
        NaCoolantChannelfvPatchScalarField
        (
            const NaCoolantChannelfvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy setting internal field reference
        NaCoolantChannelfvPatchScalarField
        (
            const NaCoolantChannelfvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new NaCoolantChannelfvPatchScalarField(*this, iF)
            );
        }


    // Member functions
    
        //- Return true if this patch field fixes a value.
        //  Needed to check if a level has to be specified while solving
        //  Poissons equations.
        virtual bool fixesValue() const
        {
            return true;
        }
        
        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Return the matrix diagonal coefficients corresponding to the evaluation of the value of this patchField
        virtual tmp<Field<scalar> > valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix source coefficients corresponding to the evaluation of the value of this patchField,
        virtual tmp<Field<scalar> > valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        //- Return the matrix diagonal coefficients corresponding to the evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientInternalCoeffs() const;

        //- Return the matrix source coefficients corresponding to the evaluation of the gradient of this patchField
        virtual tmp<Field<scalar> > gradientBoundaryCoeffs() const;
        
        //- Write to an output stream
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
