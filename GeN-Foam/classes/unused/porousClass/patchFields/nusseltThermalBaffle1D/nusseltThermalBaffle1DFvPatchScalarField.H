/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::compressible::nusseltThermalBaffle1DFvPatchScalarField


Description
    This BC simulates transient heat flux transfer through a thin baffle
    and specified heat transfer coefficients between fluid cell centers and 
    patch face in the form h = Nu*k_fluid/Dh with Nu = A + B*(Re^C)*Pr(^D). 
    The thermal inertial of the baffle is neglected. To ensure proper convergence
    within time-step, multiple PIMPLE iterations are required, yet 5 suffice for most
    cases of interest. The solid properties are specify as dictionary. Baffle and solid 
    properties need to be specified on the master side of the baffle.
    
    By Stefan Radman, stefan.radman@epfl.ch

Usage
    Example of the boundary condition specification using constant
    solid thermo :

    \verbatim
    <masterPatchName>
    {
        type            compressible::nusseltThermalBaffle1D<hConstSolidThermoPhysics>;
        value           uniform 668; // Starting T of the patch
        QInfo           0; // 1 if you want extra info print on screen regarding heat flux balance
        fh              1; // increase the Nusselt by an artificial factor fh (e.g., to assess impact of correlation uncertainty)
        thickness       uniform 0.0045; // baffle thickness
        specie
        {
            molWeight       55.845; // no effect
        }
        equationOfState
        {
            rho             8000; // no effect
        }
        thermodynamics
        {
            Cp              600; // no effect
            Hf              0; //no effect
        }
        transport
        {
            kappa           20; // thermal conductivity of the baffle
        }
        singleDh        0; // if 1, a single hydraulic diameter under Dh (to be provided, not Dh1 or Dh2) is used by both master and slave
        cellZones1      1(fuel); // list of cell zones that, if adjacent to this patch, have a hydraulic diameter of Dh1
        Dh1             0.004354644; // hydraulic diameter of those patch sides that neighbour any cellZone in cellZones1 
        cellZones2      1(gap); // list of cell zones that, if adjacent to this patch, have a hydraulic diameter of Dh2
        Dh2             0.009098158; // hydraulic diameter of those patch sides that neighbour any cellZone in cellZones2
        lowRe           1000; // laminar-to-transition Reynolds
        highRe          2300; // transition-to-turbulent Reynolds
        lNuConst        4.36; // A in Nu = A + B*(Re^C)*(Pr^D), laminar case
        lRePrCoeff      0; // B in Nu = A + B*(Re^C)*(Pr^D), laminar case
        lExpRe          0; // C in Nu = A + B*(Re^C)*(Pr^D), laminar case
        lExpPr          0; // D in Nu = A + B*(Re^C)*(Pr^D), laminar case
        tNuConst        4.6; // A in Nu = A + B*(Re^C)*(Pr^D), turbulent case
        tRePrCoeff      0.018; // B in Nu = A + B*(Re^C)*(Pr^D), turbulent case
        tExpRe          0.8; // C in Nu = A + B*(Re^C)*(Pr^D), turbulent case
        tExpPr          0.8; // D in Nu = A + B*(Re^C)*(Pr^D), turbulent case
        UName           realU; // name of physical velocity field in the solver
    }
    <slavePatchName> // The slave read almost all of the info from the master
    {
        type            compressible::nusseltThermalBaffle1D<hConstSolidThermoPhysics>;
        value           uniform 668;
        QInfo           0;
        fh              1;
    }
    \endverbatim

SourceFiles
    nusseltThermalBaffle1DFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nusseltThermalBaffle1DFvPatchScalarField_H
#define nusseltThermalBaffle1DFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "autoPtr.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

/*---------------------------------------------------------------------------*\
            Class nusseltThermalBaffle1DFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

template<class solidType>
class nusseltThermalBaffle1DFvPatchScalarField
:
    public mappedPatchBase,
    public mixedFvPatchScalarField
{
    // Private data

        //-//-//

        bool singleDh_;

        // Print additional info on heat transfer through patch
        bool QInfo_;

        // h multiplier, just for testing;
        scalar fh_;

        //- Name of the velocity field for the computation of the Reynolds
        word UName_;

        //- Names of the cellZones at the two edges of the baffle
        const List<word> cellZones1_;
        const List<word> cellZones2_;

        //- Hydraulic diameters for the computation of the heat transfer coefficient
        const scalar Dh1_;
        const scalar Dh2_;

        //- Actual patch hydraulic diameter. Set by setProperties().
        mutable scalar Dh_;

        // Reynolds threshold number for laminar-to-transition
        // and transition-to-turbulent regions
        mutable scalar lowRe_;
        mutable scalar highRe_;

        //- Coefficients in Nu = NuConst + RePrCoeff*Re^(expRe)*Pr^(expPr)
            //- Laminar
        mutable scalar lNuConst_;
        mutable scalar lRePrCoeff_;
        mutable scalar lExpRe_;
        mutable scalar lExpPr_;

            //- Turbulent
        mutable scalar tNuConst_;
        mutable scalar tRePrCoeff_;
        mutable scalar tExpRe_;
        mutable scalar tExpPr_;

        //- Real baffle temperature at interface with fluid
        mutable scalarField TpReal_;

        //- Heat transfer coefficient between the patch center and first cell center
        mutable scalarField h_;

        //- Internal temperature field of the first cells next to patch
        mutable scalarField Ti_;

        //- Field btw 0 and 1 to select whether I am laminar or turbulent
        mutable scalarField fRe_;

        mutable scalarField Nu_;

        //-//-//

        //- Name of the temperature field
        word TName_;

        //- Baffle is activated
        bool baffleActivated_;

        //- Baffle thickness [m]
        mutable scalarField thickness_;

        //- Solid dictionary
        dictionary solidDict_;

        //- Solid thermo
        mutable autoPtr<solidType> solidPtr_;


    // Private members

        //-//-//

        //- Set the correct hydraulic diameter of the patch and the 
        //- Nusselt correlation parameters
        void setProperties();

        //- Compute the heat transfer coefficient for this side of the patch()
        void computeH();

        //- Compute the fraction of turbulent flow (0 = lam, 1 = turb, in between = transition)
        void computeFRe();
        
        //- Return temperature field of neighbour;
        tmp<scalarField> nbrTpReal() const;

        //- Return heat transfer coefficient for the other side of the patch;
        tmp<scalarField> nbrH() const;

        //- Return internalPatchField temperature of neighbour;
        tmp<scalarField> nbrTi() const;

        //-//-//

        //- Return const solid thermo
        const solidType& solid() const;

        //- Return thickness from master
        tmp<scalarField> baffleThickness() const;

        //- Is Owner
        bool owner() const;


public:

    //- Runtime type information
    TypeName("compressible::nusseltThermalBaffle1D");


    // Constructors

        //- Construct from patch and internal field
        nusseltThermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nusseltThermalBaffle1DFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  nusseltThermalBaffle1DFvPatchScalarField onto a new patch
        nusseltThermalBaffle1DFvPatchScalarField
        (
            const nusseltThermalBaffle1DFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nusseltThermalBaffle1DFvPatchScalarField
        (
            const nusseltThermalBaffle1DFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nusseltThermalBaffle1DFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nusseltThermalBaffle1DFvPatchScalarField
        (
            const nusseltThermalBaffle1DFvPatchScalarField&,
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
                new nusseltThermalBaffle1DFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "nusseltThermalBaffle1DFvPatchScalarField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
