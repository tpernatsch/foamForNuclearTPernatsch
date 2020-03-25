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

    //- Precursor "power" for each delayed group in W to be compatible with
    //  point kinetics being solved for power rather than neutron density
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
    scalar coeffTClad_;
    scalar coeffTCool_;
    scalar coeffRhoCool_;
    scalar coeffTStruct_;

    //- Reference integral values for feedbacks
    scalar TFuelRef_;
    scalar TCladRef_;
    scalar TCoolRef_;
    scalar rhoCoolRef_;
    scalar TStructRef_;

    //- Ptrs to coupling fields (i.e. defined on their original mesh)
    const volScalarField* TFuelOrig_;
    const volScalarField* TCladOrig_;
    const volScalarField* TCoolOrig_;
    const volScalarField* rhoCoolOrig_;
    const volScalarField* TStructOrig_;
    const volScalarField* volFuelPowerOrig_;

    //- These have no use for now, will be needed by the liquidFuel version
    const volVectorField* UOrig_;
    const volScalarField* alphaOrig_;
    const volScalarField* alphatOrig_;
    const volScalarField* muOrig_;

    //- Mapped fields
    volScalarField TFuel_;
    volScalarField TClad_;
    volScalarField TCool_;
    volScalarField rhoCool_;
    volScalarField TStruct_;

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
    volScalarField coolFeedbackCellField_;

    //- For TStructures
    volScalarField structFeedbackCellField_;

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

        //- Return fluxes
        virtual const PtrList<volScalarField> fluxes() const
        {
            return fluxes_;
        };

        //- Return precursors
        virtual const PtrList<volScalarField> precursors() const
        {
            return precursors_;
        };

        //- Return one group flux
        virtual const volScalarField& oneGroupFlux() const
        {
            return oneGroupFlux_;
        }

        //- Return power
        virtual scalar power() const
        {
            return power_;
        };

        virtual void getCouplingFieldRefs
        (
            const objectRegistry& src,
            const meshToMesh& neutroToFluid
        );   

        virtual void interpolateCouplingFields
        (
            const meshToMesh& neutroToFluid
        );

        //- Correct/update the properties
        virtual void correct
        (
            scalar& residual,
            label couplingIter  
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
