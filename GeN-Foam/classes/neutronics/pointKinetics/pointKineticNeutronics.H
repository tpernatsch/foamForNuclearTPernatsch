/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::neutronics::pointKineticNeutronics

Description
    Run-time selectable point-kinetic based neutronics solver

SourceFiles
    pointKineticNeutronics.C

Authors
    Stefan Radman <stefanradman92@gmail.com; stefan.radman@epfl.ch;>
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

\*---------------------------------------------------------------------------*/

#ifndef pointKineticNeutronics_H
#define pointKineticNeutronics_H

#include "neutronics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class pointKineticNeutronics Declaration
\*---------------------------------------------------------------------------*/

class pointKineticNeutronics
:
    public neutronics
{
protected:

    //-
    IOdictionary nuclearData_;

    //-
    scalar power_;

    //- Due to whatever the user might please
    scalar externalReactivity_;

    //-
    scalar totalReactivity_;

    //-
    scalar promptGenerationTime_;
    
    //- Total effective delayed neutron fraction
    scalar beta_;

    //- Effective delayed fractions for each delayed group
    scalarList betas_;

    //- Decay constants for each delayed group
    scalarList lambdas_;

    //- Precursor "power" for each delayed group
    //  in W to be compatible with point kinetics
    //  solved for power rather than neutron density
    scalarList precursorPowers_;

    //- Number of delayed groups
    label delayedGroups_;

    //- Old time power, precursorPowers to be used to calculate derivatives,
    //  as well as timeIndex to check if time step has changed (can't think
    //  of a better way to get something equivalent to the finalInnerIter flag)
    label timeIndex_;
    scalar powerOld_;
    scalarList precursorPowersOld_;

    //- Feedback coefficients
    scalar coeffTFuel_;
    scalar coeffTCladding_;
    scalar coeffTCoolant_;
    scalar coeffRhoCoolant_;

    //- Adaptations for the FFTF, all radial/axial structure expansion
    //  feedbacks indirectly modelled via the the corresponding structure
    //  temperature change that caused them
    scalar coeffTStructures_;

    //- Reference integral values for feedbacks
    scalar TFuelRef_;
    scalar TCladdingRef_;
    scalar TCoolantRef_;
    scalar rhoCoolantRef_;
    scalar TStructuresRef_;

    //- If the Ref values are present in nuclearData, read those and do not
    //  compute them in the first correct. I need this because I cannot get
    //  fields at construction to compute Ref values if they are not provided.
    //  Thus, if they are provided, I need to find a way to tell the correct
    //  not to recompute them. This is the only way for now (very ugly...)
    bool resetTFuelRef_;
    bool resetTCladdingRef_;
    bool resetTCoolantRef_;
    bool resetRhoCoolantRef_;
    bool resetTStructuresRef_;

    //- Fields for feedbacks
    volScalarField TFuel_;
    volScalarField TCladding_;
    volScalarField TCoolant_;
    volScalarField rhoCoolant_;
    volScalarField TStructures_;

    //- Well, these need to be passed via getFields but, ooops, getField
    //  does not support those fields. This needs to be changed in all the
    //  neutronic classes, not doing it now
    //volScalarField axialFuelExpansion_;
    //volScalarField coreRadialExpansion_;
    //volScalarField claddingRadialExpansion_;

    //-
    volScalarField oneGroupFlux_;

    //- Number of energy groups (for scaling fluxes of a possibly existing
    //  previous simulation)
    label energyGroups_;

    //- Group-wise fluxes
    PtrList<volScalarField> fluxes_;

    //- Delayed-group-wise precursor concentrations
    PtrList<volScalarField> precursors_;

    //- Delayed-group-wise precursor power fields, to be used for MSRs
    PtrList<volScalarField> liquidFuelPrecursorPowers_;

    //- If true, fuel temperature feedback (Doppler) computed as 
    //  deltaReactivity = -coeffFastDoppler*ln(TFuel/TFuelRef) with 
    //  coeffFastDoppler being the Doppler constant. If false, Doppler feedback
    //  computed as deltaReactivity = coeffTFuel*(TFuel-TFuelRef)
    bool fastNeutrons_;

    //- Doppler constant if doing fastNeutrons
    scalar coeffFastDoppler_;

    //- Fields to mark which cells are part of X and which aren't, used for
    //  computing the 0-D feedback values (TFuel, TStructures, etc...).
    //  Indeed, the spatial extent of these marker fields should reflect
    //  the spatial extent of the domain that was perturbed in neutronics
    //  calculations to compute feedback coefficients
    
    //- For TFuel and TCladding
    volScalarField fuelFeedbackCellField_;

    //- For rhoCoolant and TCoolant
    volScalarField coolantFeedbackCellField_;

    //- For TStructures
    volScalarField structuresFeedbackCellField_;

    //-----------------------------------------------------------------------//
    //- For the FFTF only, GEM and driveline --------------------------------//
    //-----------------------------------------------------------------------//
    
    //- Coolant flow
    surfaceScalarField phi_;

    //- List of face labels over which to compute integral flow values
    labelList flowFaces_;

    //- Reference value of the integral flow value
    scalar flowRef_;
    
    //- Core flow as a fraction of the reference flow
    scalar flowFraction_; 

    //- Reference flow fraction at the start of the simulation
    scalar flowFractionRef_;

    //- An optional offset to the sodium height for fine adjustments
    //  without having to change all entries in GEMReactivityMap
    scalar sodiumLevelRef_;

    //- Pairs of GEM sodium level (in m wrt a reference) vs reactivity (in
    //  absolute units)
    List<Pair<scalar>> GEMReactivityMap_;

    //- 
    volScalarField drivelineFeedbackCellField_;

    //- Reference driveline temperature
    scalar TDrivelineRef_;

    //- Absolute expansion coeff of the driveline in m/K
    scalar drivelineExpansionCoeff_;

    //- Pairs of CR insertion value (in m wrt a reference) vs reactivity (in
    //  absolute)
    List<Pair<scalar>> CRReactivityMap_;

    //- For ease of printing
    HashTable<scalar, word, word::hash> reactivityValues_;

private:

    //- Disallow default bitwise copy construct
    pointKineticNeutronics(const pointKineticNeutronics&);

    //- Disallow default bitwise assignment
    void operator=(const pointKineticNeutronics&);

public:

    //- Runtime type information
        TypeName("pointKinetics");

    //- Constructors

        //- Construct from mesh, read data from IOdictionaries
        pointKineticNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~pointKineticNeutronics(){}

    //- Member Functions

        //- Solve pointKinetics
        virtual void correct
        (
            const label couplingIter, 
            scalar& residual, 
            const bool& liquidFuel
        );  

        //- Take fields to compute 0-D reactivity feedbacks
        virtual void getFields
        (
            const volVectorFieldTable& vvfTable,
            const volScalarFieldTable& vsfTable, 
            const meshToMesh& neutroToFluid
        );  

        //- "Override" getFieldsLiquidFuel
        virtual void getFieldsLiquidFuel
        (
            const volVectorFieldTable& vvfTable,
            const volScalarFieldTable& vsfTable, 
            const meshToMesh& neutroToFluid
        );

        //- Return fluxes
        virtual inline const PtrList<volScalarField> fluxes() const
        {
            return fluxes_;
        };

        //- Return precursors
        virtual inline const PtrList<volScalarField> precursors() const
        {
            return precursors_;
        };

        //- Return one group flux
        virtual inline const volScalarField& Flux() const
        {
            return oneGroupFlux_;
        }

        //- Return power
        virtual scalar power() const
        {
            return power_;
        };

        //- For the FFTF only

            //- Gas Expansion Modules
            scalar GEMReactivity(); 

            //- Driveline thermal expansion
            scalar drivelineReactivity();

            //- Print individual reactivity contributions
            virtual void printInfo() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
