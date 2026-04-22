/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
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
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::coolantChannelModel

Description
    Coolant channel state model used by the coolantChannel boundary condition.

    The class stores channel properties, coolant boundary inputs and derived
    coolant thermodynamic state. Patch-specific geometry quantities, such as
    face heights and axial face locations, are supplied by the boundary
    condition during the update.

SourceFiles
    coolantChannelModel.C

\*---------------------------------------------------------------------------*/

#ifndef coolantChannelModel_H
#define coolantChannelModel_H

#include "scalar.H"
#include "scalarField.H"
#include "dictionary.H"
#include "autoPtr.H"

#ifdef OPENFOAMESI
    #include "Function1.H"
#endif

#include "Table.H"
#include "InterpolateTables.H"
#include "scalarFieldFieldINew.H"
#include "PtrListScalarFieldFieldINew.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class coolantChannelModel Declaration
\*---------------------------------------------------------------------------*/

class coolantChannelModel
{
private:

    // Private data

        //- Coolant flow area
        scalar A_;

        //- Flow direction along pin axis (+1 or -1)
        scalar flowDirection_;

        //- Hydraulic diameter
        scalar DH_;

        //- Coolant pressure
        scalarField coolantPressure_;

#ifdef OPENFOAMFOUNDATION
        typedef Function1s::Table<scalar> scalarTable;

        //- Time-dependent coolant pressure
        autoPtr<scalarTable> coolantPressureList_;
#elif OPENFOAMESI
        typedef Foam::Function1Types::Table<scalar> scalarTable;

        //- Time-dependent coolant pressure
        autoPtr<scalarTable> coolantPressureList_;
#endif

        //- Mass flow rate
        scalarField massFlowRate_;

#ifdef OPENFOAMFOUNDATION
        //- Time-dependent mass flow rate
        autoPtr<scalarTable> massFlowRateList_;
#elif OPENFOAMESI
        //- Time-dependent mass flow rate
        autoPtr<scalarTable> massFlowRateList_;
#endif

        //- Inlet coolant temperature
        scalar inletTemperature_;

#ifdef OPENFOAMFOUNDATION
        //- Time-dependent inlet temperature
        autoPtr<scalarTable> inletTemperatureList_;
#elif OPENFOAMESI
        //- Time-dependent inlet temperature
        autoPtr<scalarTable> inletTemperatureList_;
#endif

        //- If true, compute bulk state through enthalpy integration
        bool enthalpyModel_;

        //- If true, inlet temperature is required by another model
        bool inletTemperatureRequired_;

        //- If true, all mass flow rate values must be positive
        bool positiveMassFlowRequired_;

        //- Axial coordinates for prescribed temperature profile
        scalarField zValues_;

        //- Time coordinates for prescribed temperature profile
        scalarField timeValues_;

        //- Temperature profile data
        FieldField<Field, scalar> T0ProfileData_;

        //- Axial interpolation method
        interpolateTableBase::interpolationMethod zMethod_;

        //- Time interpolation method
        interpolateTableBase::interpolationMethod tMethod_;

        //- Time interpolation table for axial temperature profiles
        autoPtr<scalarFieldInterpolateTable> T0Table_;

        //- Stored axial profile dictionary
        dictionary axialProfileDict_;

        //- Bulk coolant temperature at patch faces
        scalarField T0_;

        //- Steam quality based on enthalpy
        scalarField Q_;

        //- Void fraction estimate
        scalarField X_;

        //- Reynolds number
        scalarField Re_;

        //- Prandtl number
        scalarField Pr_;

        //- Film temperature
        scalarField Tfilm_;

        //- Coolant thermal conductivity
        scalarField coolantk_;

        //- Coolant specific heat
        scalarField coolantCp_;

        //- Coolant enthalpy
        scalarField coolantEnthalpy_;

        //- Mixture density at each face [kg/m³] — transient state
        scalarField rho_;

        //- Mixture density at previous time step [kg/m³]
        scalarField rhoOld_;

        //- Specific enthalpy at previous time step [J/kg]
        scalarField hOld_;

        //- Mixture velocity at each face [m/s]
        scalarField velocity_;

        //- Maximum Picard iterations per cell per updateBulkState call
        label maxPicardIter_;

        //- Enthalpy convergence tolerance for Picard iteration [J/kg]
        scalar picardTol_;

        //- Time tracker for freezing old-time state in HEM (initialized to -GREAT)
        scalar timeCheckHEM_;

        //- Inlet mass flow rate [kg/s] — preserved across HEM overwrites of massFlowRate_
        scalar mdotInlet_;


    // Private Member Functions

        //- Initialize prescribed axial temperature profile
        void initializeAxialProfile(const dictionary& dict);

        //- Rebuild the prescribed axial temperature profile table
        void initializeAxialProfileTable();

        //- Validate channel properties
        void validate() const;


public:

    // Constructors

        //- Construct from number of faces
        coolantChannelModel(const label nFaces);

        //- Construct from dictionary and number of faces
        coolantChannelModel
        (
            const dictionary& dict,
            const label nFaces,
            const bool inletTemperatureRequired = false,
            const bool enthalpyModelDefault = true,
            const bool positiveMassFlowRequired = true
        );

        //- Copy constructor
        coolantChannelModel(const coolantChannelModel& model);


    // Member Functions

        //- Update time-dependent boundary inputs
        void updateBoundaryInputs(const scalar userTime);

        //- Compute bulk coolant temperature and enthalpy along the channel.
        //  deltaT is the OpenFOAM time step size [s]; defaults to GREAT
        //  which gives the steady-state limit (backward-compatible).
        void updateBulkState
        (
            const scalar userTime,
            const scalarField& zetas,
            const scalarField& radii,
            const scalarField& faceHeights,
            const scalarField& heatFlux,
            const scalar deltaT = GREAT
        );

        //- Compute derived coolant thermodynamic state
        void updateDerivedState(const scalarField& Twall);

        //- Write model entries and coolant state
        void write(Ostream& os) const;


    // Access

        scalar A() const
        {
            return A_;
        }

        scalar flowArea() const
        {
            return A_;
        }

        scalar flowDirection() const
        {
            return flowDirection_;
        }

        scalar DH() const
        {
            return DH_;
        }

        scalar hydraulicDiameter() const
        {
            return DH_;
        }

        scalar inletTemperature() const
        {
            return inletTemperature_;
        }

        bool enthalpyModel() const
        {
            return enthalpyModel_;
        }

        const scalarField& T0() const
        {
            return T0_;
        }

        const scalarField& coolantPressure() const
        {
            return coolantPressure_;
        }

        const scalarField& pressure() const
        {
            return coolantPressure_;
        }

        const scalarField& massFlowRate() const
        {
            return massFlowRate_;
        }

        const scalarField& coolantEnthalpy() const
        {
            return coolantEnthalpy_;
        }

        const scalarField& enthalpy() const
        {
            return coolantEnthalpy_;
        }

        const scalarField& Q() const
        {
            return Q_;
        }

        const scalarField& quality() const
        {
            return Q_;
        }

        const scalarField& X() const
        {
            return X_;
        }

        const scalarField& voidFraction() const
        {
            return X_;
        }

        const scalarField& Re() const
        {
            return Re_;
        }

        const scalarField& Pr() const
        {
            return Pr_;
        }

        const scalarField& Tfilm() const
        {
            return Tfilm_;
        }

        const scalarField& coolantk() const
        {
            return coolantk_;
        }

        const scalarField& coolantCp() const
        {
            return coolantCp_;
        }

        const scalarField& velocity() const
        {
            return velocity_;
        }

        const scalarField& rho() const
        {
            return rho_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
