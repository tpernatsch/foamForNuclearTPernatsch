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
    Foam::coolantChannelfvPatchScalarField

Description
    Coolant channel boundary condition for thermal coupling between a solid wall
    and a 1D coolant stream represented through a patch.

    The model computes:
    - wall heat flux from the solid side,
    - axial bulk coolant temperature/enthalpy,
    - derived coolant state quantities,
    - heat transfer coefficients using selected correlations for different
      regimes,
    - mixed boundary coefficients for the wall temperature.

Usage
    Example:
    \verbatim
        cladOuter
        {
            type                    coolantChannel;
            value                   uniform 560;

            kappa                   k;
            hydraulicDiameter       0.01;
            flowArea                1e-4;
            flowDirection           1.0;
            enthalpyModel           true;
            relax                   0.15;

            correlations
            {
                singlePhase         Gnielinski;
                onb                 Basu;
                chf                 EPRI;
                subcooledBoiling    Chen;
                saturatedBoiling    Chen;
                transitionBoiling   McDonoughMilichKing;
                filmBoiling         Frederking;
                leidenfrost         GroeneveldStewartCorrected;
            }

            inletTemperature        560;
            coolantPressure         uniform 1.55e7;
            massFlowRate            uniform 0.25;
        }
    \endverbatim

    Time-dependent inlet data can also be provided with:
    - inletTemperatureList
    - coolantPressureList
    - massFlowRateList

    If enthalpyModel is false, an axialProfileDict can be provided instead.

SourceFiles
    coolantChannelfvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef coolantChannelfvPatchScalarField_H
#define coolantChannelfvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

#include "sliceMapper.H"

#include "coolantChannelModel/coolantChannelModel.H"
#include "singlePhaseHTCModel.H"
#include "onbModel.H"
#include "subcooledBoilingHTCModel.H"
#include "saturatedBoilingHTCModel.H"
#include "filmBoilingHTCModel.H"
#include "chfModel.H"
#include "leidenfrostModel.H"
#include "transitionBoilingHTCModel.H"

namespace Foam
{

class coolantChannelfvPatchScalarField
:
    public fvPatchField<scalar>
{
private:

    // Private data

        //- Name of thermal conductivity field in the solid region
        word kappaName_;

        //- Correlation selection dictionary
        dictionary correlations_;

        // Heat transfer regime models
        singlePhaseHTCModel    singlePhase_;
        onbModel               onb_;
        subcooledBoilingHTCModel subcooledBoiling_;
        saturatedBoilingHTCModel saturatedBoiling_;
        filmBoilingHTCModel    filmBoiling_;
        chfModel               chf_;
        leidenfrostModel       leidenfrost_;
        transitionBoilingHTCModel transitionBoiling_;

        //- Coolant channel properties and state
        coolantChannelModel coolant_;

        //- Heat transfer coefficient at patch faces
        scalarField h_;

        //- Mixed BC coefficient
        scalarField alpha_;

        //- Relaxation factor for HTC update
        scalar relax_;

        //- Auxiliary time marker used to control updates
        scalar timeCheck_;

        //- Wall heat flux
        scalarField heatFlux_;

        //- Face heights used for axial integration of the coolant energy balance
        scalarField faceHeights_;


    // Private Member Functions

        //- Compute patch face heights projected along the pin direction
        scalarField calculateFaceHeights() const;

        //- Validate geometry and required user inputs
        void validateGeometry() const;

        //- Update face heights when needed
        void updateFaceHeightsIfNeeded();

        //- Compute solid-side conductance alphaP
        scalarField computeSolidSideConductance() const;

        //- Compute wall heat flux from current wall temperature estimate
        void computeHeatFlux(const scalarField& alphaP);

        //- Solve pre-CHF wall temperature for a target heat flux.
        //  Self-consistent: finds T_w such that q_SB(T_w) = qTarget.
        scalar solvePreCHFTemperature
        (
            const label faceI,
            const scalar qTarget
        ) const;

        //- Solve film-boiling wall temperature for a target heat flux.
        //  Self-consistent: finds T_w such that q_FB(T_w) = qTarget.
        scalar solveFilmBoilingTemperature
        (
            const label faceI,
            const scalar qTarget
        ) const;


public:

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
            const bool valueRequired = true
        );

        //- Construct by mapping given field onto a new patch
        coolantChannelfvPatchScalarField
        (
            const coolantChannelfvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired = true
        );

        //- Construct as copy setting internal field reference
        coolantChannelfvPatchScalarField
        (
            const coolantChannelfvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new coolantChannelfvPatchScalarField(*this)
            );
        }

        //- Clone setting internal field reference
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


    // Member Functions

        virtual bool fixesValue() const
        {
            return true;
        }

        virtual void updateCoeffs();

        virtual void evaluate
        (
            const Pstream::commsTypes commsType = Pstream::commsTypes::blocking
        );

        virtual tmp<Field<scalar>> valueInternalCoeffs
        (
            const tmp<scalarField>&
        ) const;

        virtual tmp<Field<scalar>> valueBoundaryCoeffs
        (
            const tmp<scalarField>&
        ) const;

        virtual tmp<Field<scalar>> gradientInternalCoeffs() const;

        virtual tmp<Field<scalar>> gradientBoundaryCoeffs() const;

        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
