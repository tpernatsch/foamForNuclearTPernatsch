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
    Foam::coolantChannelfvPatchScalarField

Description
    Coolant subchannel module to provide thermal boundary condition.
    
    The module can be used in two modes. Either the time-dependent axial bulk
    coolant temperature profile is provided or the enthalpy model must be
    enabled and the constant or time-dependent inlet coolant temperature, 
    mass flow rate and pressure data must be supplied. The former mode can 
    model natural and forced flow as well, while the latter, by principle,
    can only be used in the case of forced coolant flow.

    The enthalpy model uses the closed channel assumption. Only the mass 
    and energy conservation equations are effectively solved. Pressure drop
    and other inlet/outlet effects are neglected. Futhermore, the model uses
    the homogeneous equilibrium model (HEM; i.e. the slip ratio is unity).
    
    At its current state, the coolantChannel module can only be used for 
    light-water as coolant. The material properties are obtained from the
    CoolProp IF97 C++ functions of (c) Ian Bell and the CoolProp team.
    The original implementation rely on the IF97 formulation by IAPWS. 
    The actually used material library was further extended with additional 
    material property functions.
    
    **Verification**:

    The module was compared against the U.S. NRC system code TRACE for
    various simplified 1D rod cases. The investigated pins included 
    traditional and novel PWR and BWR designs with forced coolant flow. 
    The test scenarios were designed to cover the whole range of the
    boiling curve with the exception of transition boiling and single phase
    vapor flow. The module is considered to be successfully verified until
    departure from nucleate boiling (DNB) occurs. Post-DNB wall temperatures
    (especially for PWR cases) showed larger than expected deviations from
    the reference TRACE results.
________________________________________________________________________________
\par Options
    The `coolantChannel` fvPatchField can be selected in the patch subdictionary
    inside the `boundaryField` subdictionary of the temperature field.

    In case CHF model `lookuptableGroeneveld` is selected, the user must provide
    CHF table as IOdictionary in `constant/chfTableDict`.

    The dictionary `constant/chfTableDict` should look like:
    \verbatim
        pressureValues 
        (
            100.0   300.0  ... 21000.0
        );

        massFluxValues 
        (
            0.0    50.0 ... 8000.0
        );

        qualityValues
        (
            -0.50 -0.40 ... 1.0
        ); 

        data 
        (
            // For p = 100.0                                                                                         
            (                                                                                           
                (8111      7252    ...      0)       // For G = 0.0
                (44338     36803   ...     0)        // For G = 8000.0
            )                                                                                           
            // For p = 300.0                                                                                            
            (
                ...                                                                              
            )                                                                                           
            ...                                                                                           
            // For p = 21000.0                                                                                          
            (
                ...                                                                                           
            )                                                                                           
        );
    \endverbatim
________________________________________________________________________________
Usage
    **Case 1 :** Enthalpy model disabled, bulk coolant temperature provided.
    \verbatim
        cladOuter
        {
            type                coolantChannel;
            value               uniform 300.0;
            hydraulicDiameter   0.01;
            enthalpyModel       0;
            axialProfileDict
            {
                timePoints      ( 0.0 3600.0 64000.0 );
                axialLocations  ( 0 0.5 1.0 1.5 );
                T0Data          ( ( 293.0 293.0 293.0 293.0 ) 
                                ( 440.0 450.0 460.0 460.0 ) 
                                ( 475.0 490.0 500.0 510.0 ) );
            }
            relaxHTC            0.7;
            flowArea            0.0001;
            correlations
            {
                ONB                 Basu;
                CHF                 EPRI;
                Leidenfrost         GroeneveldStewartCorrected;
                SinglePhase         Gnielinski;
                SubcooledBoiling    Rohsenow;
                SaturatedBoiling    Rohsenow;
                TransitionBoiling   McDonoughMilichKing;
                FilmBoiling         GroeneveldTubular;
            }
            coolantPressure         uniform 1e7;
            massFlowRate            uniform 0.25;
        }

    \endverbatim

    **Case 2 :** Enthalpy model enabled, inlet temperature provided.
    \verbatim
        cladOuter
        {
            type                coolantChannel;
            value               uniform 300.0;
            hydraulicDiameter   0.01;
            enthalpyModel       1;
            relaxHTC            0.7;
            flowArea            0.0001;
            flowDirection       1.0;
            correlations
            {
                ONB                 Basu;
                CHF                 lookuptableGroeneveld;
                Leidenfrost         GroeneveldStewartCorrected;
                SinglePhase         Gnielinski;
                SubcooledBoiling    Rohsenow;
                SaturatedBoiling    Rohsenow;
                TransitionBoiling   McDonoughMilichKing;
                FilmBoiling         GroeneveldTubular;
            }
            inletTemperature        300;
            coolantPressure         uniform 101325;
            massFlowRate            uniform 0.25;
        }
    \endverbatim

SourceFiles
    coolantChannelfvPatchScalarField.C

\mainauthor
    A. Sz. Vanyi - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    A. Scolaro, E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\date 
    July 2024

\*---------------------------------------------------------------------------*/

#ifndef coolantChannelfvPatchScalarField_H
#define coolantChannelfvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "axialProfile.H"
#include "InterpolateTables.H"
#ifdef OPENFOAMESI
    #include "Function1.H"
#endif     
#include "Table.H"
#include "sliceMapper.H"
#include "PtrListScalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class coolantChannelfvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/
class coolantChannelfvPatchScalarField
:
    public fvPatchField<scalar>
{
    // Private data
        //- Name of thermal conductivity field
        word kappaName_;

        //- Dictionary with user-defined correlations
        dictionary correlations_;

        //- List of HTC correlations for the axial slices
        List<word> correlationRegimeList_;
        
        //- T0 patch profile @ current timestep - working variable
        scalarField T0_;

        //- h patch profile @ current timestep - for write()
        scalarField h_;
        
        //- Boundary heat transfer factor - working variable
        scalarField alpha_;

        //- Subchannel flow area
        scalar A_;

         //- Direction of the forced flow
        scalar flowDirection_;

        //- Relaxation factor of the HTC
        scalar relaxHTC_;

        //- Film boiling correction factor for RIA transients
        scalar filmBoilingRIACorrection_;
        
        //- Transition boiling correction factor for RIA transients
        scalar transitionBoilingRIACorrection_;

        //- Temperature increment to define the rewetting temperature for RIA transients 
        scalar leidenfrostdTBessiron_;
        
        //- Water layer thickness to be vaporized for RIA transients
        scalar thicknessBessironRIA_;

        //- Auxiliary variable for controlling vaporization in case of RIA transients
        scalar timeCheck_;

        //- Hydraulic diameter
        scalar DH_;

        //- Critical heat flux
        scalarField qCHF_;

        //- Vaporized thickness for RIA transients
        scalarField vaporizedThickness_;

        //- Wall heat flux
        scalarField heatFlux_;

        //- Steam quality
        scalarField Q_;

        //- Void fraction
        scalarField X_;

        //- Reynolds number
        scalarField Re_;

        //- Prandtl number
        scalarField Pr_;

        //- Wall temperature
        scalarField Twall_;

        //- Film temperature
        scalarField Tfilm_;

        //- Coolant thermal conductivity
        scalarField coolantk_;

        //- Coolant enthalpy
        scalarField coolantEnthalpy_;

        //- Coolant pressure
        scalarField coolantPressure_;

#ifdef OPENFOAMFOUNDATION
        //- Coolant pressure time-dependent list
        autoPtr<Function1s::Table<scalar>> coolantPressureList_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> scalarTable;
        //- Coolant pressure time-dependent list
        autoPtr<scalarTable> coolantPressureList_;
#endif

        //- Coolant mass flow rate
        scalarField massFlowRate_;

#ifdef OPENFOAMFOUNDATION
        //- Coolant mass flow rate time-dependent list
        autoPtr<Function1s::Table<scalar>> massFlowRateList_;
#elif OPENFOAMESI
        //- Coolant mass flow rate time-dependent list
        autoPtr<scalarTable> massFlowRateList_;
#endif     

        //- Inlet temperature
        scalar inletTemperature_;

#ifdef OPENFOAMFOUNDATION
        //- Inlet temperature time-dependent list
        autoPtr<Function1s::Table<scalar>> inletTemperatureList_;
#elif OPENFOAMESI
        //- Inlet temperature time-dependent list
        autoPtr<scalarTable> inletTemperatureList_;
#endif 

        //- Flag for enthalpy model
        bool enthalpyModel_;
        
        //- Tabulated CHF data
        PtrList<FieldField<Field, scalar>> chfData_;
        
        //- Optional P/D to CHF correction
        scalar chfPitchToDiameter_;
        
        //- Pressure list for CHF tables
        scalarField chfPressureValues_;
        
        //- Mass flux list for CHF tables
        scalarField chfMassFluxValues_;
        
        //- Quality list for CHF tables
        scalarField chfQualityValues_;
        
        //- Interpolation method for CHF table (pressure)
        interpolateTableBase::interpolationMethod chfPressureMethod_;
        
        //- Interpolation method for CHF table (mass flux)
        interpolateTableBase::interpolationMethod chfMassFluxMethod_;
        
        //- Interpolation method for CHF table (quality)
        interpolateTableBase::interpolationMethod chfQualityMethod_;
        
        //- 3D interpolation table for pressure interpolation
        autoPtr<scalarFieldFieldInterpolateTable> chfPressureTable_;

        //- Tabulated axial locations
        scalarField zValues_;

        //- Tabulated time values
        scalarField timeValues_;

        //- Tabulated Coolant Profile data
        FieldField<Field, scalar> T0ProfileData_;

        //- Axial location interpolation method
        interpolateTableBase::interpolationMethod zMethod_;
        
        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;

        //- T0 Interpolation table
        autoPtr<scalarFieldInterpolateTable> T0Table_;

        //- Axial profile dict
        dictionary dict_;

        //- Face heights
        scalarField faceHeights_;

        //- Unit conversion factors
        const scalar PsiToPa_;
        const scalar BtuToJ_;
        const scalar FeetToMeters_;
        const scalar LbToKg_;
        const scalar HrToSec_;
        const scalar g_;

   // Private functions

        //- Calculate face heights
        virtual scalarField calculateFaceHeights();

        //- Initialize correlations dict
        virtual dictionary initCorrelationsDict(const dictionary& dict);

public:

    //- Runtime type information
    TypeName("coolantChannel");


    // Constructors

        //- Construct from patch and internal field
        coolantChannelfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        coolantChannelfvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping given coolantChannelfvPatchScalarField
        //  onto a new patch
        coolantChannelfvPatchScalarField
        (
            const coolantChannelfvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );

        //- Construct as copy setting internal field reference
        coolantChannelfvPatchScalarField
        (
            const coolantChannelfvPatchScalarField&,
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
                new coolantChannelfvPatchScalarField(*this, iF)
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

        // Mapping functions
        
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
