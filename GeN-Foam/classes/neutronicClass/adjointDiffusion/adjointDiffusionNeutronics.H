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
    Foam::neutronics::adjointDiffusionNeutronics

Description
    subclass for adjointDiffusion neutronics.

SourceFiles
    adjointDiffusionNeutronics.C

Authors
	Euqtions: Jun Shi <junshi@berkeley.edu>
    Coding: Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

\*---------------------------------------------------------------------------*/

#ifndef adjointDiffusionNeutronics_H
#define adjointDiffusionNeutronics_H

#include "neutronics.H"
#include "XS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class adjointDiffusionNeutronics Declaration
\*---------------------------------------------------------------------------*/

class adjointDiffusionNeutronics
:
    public neutronics//class inheritance
{
protected:

    // Protected data

        // Cross sections
        XS xs_;

        //Fields
        volScalarField Dalbedo_;
        PtrList<volScalarField> flux_;
        PtrList<volScalarField> fluxStar_;
        PtrList<volScalarField> prec_;
        volScalarField fluxStarAlbedo_;
        volScalarField defaultFlux_;
        volScalarField defaultPrec_;

        //Derived fields
        volScalarField oneGroupFlux_;
        volScalarField scatteringSourceExtra_;

        //mappedFields
        volVectorField U_;
        surfaceScalarField phi_;
        volScalarField porosity_;
        volScalarField Tfuel_;
        volScalarField Tclad_;
        volScalarField rhoCool_;
        volScalarField TCool_;
        volScalarField diffCoeffPrec_;

        // calc
        List<scalar> calculated_beta_eff;
        volScalarField adjointPrec_beta_;
        volScalarField adjointFlux_chiP_;
        volScalarField adjointFlux_chiD_;

private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        adjointDiffusionNeutronics(const adjointDiffusionNeutronics&);

        //- Disallow default bitwise assignment
        void operator=(const adjointDiffusionNeutronics&);


public:

    //- Runtime type information
        TypeName("adjointDiffusionNeutronics");


    // Constructors

        //- Construct from mesh, read data from IOdictionaries
        adjointDiffusionNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~adjointDiffusionNeutronics();


    // Member Functions

        // Access

            //- Return fluxes
            virtual inline const PtrList<volScalarField> fluxes() const
            {
                return flux_;
            }

            //- Return precursors
            virtual inline const PtrList<volScalarField> precursors() const
            {
                return prec_;
            }

            //- Return one group flux
            virtual inline const volScalarField& Flux() const
            {
                return oneGroupFlux_;
            }

        // Calc
        //- Return power
            virtual inline scalar power() const
            {
                scalar pTot = 0;
                forAll(flux_,energyI)
                {
                    
                    pTot += fvc::domainIntegrate(flux_[energyI] * xs_.sigmaPow()[energyI]).value();

                }
                return pTot;
            }   

            //- Correct/update the properties
            virtual void correct(const label couplingIter, scalar& residual, const bool& liquidFuel);  

            //- Take from another mesh the fields for XS parametrization
            virtual void getFields(
            const volScalarField& TfuelOrig, 
            const volScalarField& TcladOrig, 
            const volScalarField& rhoCoolOrig, 
            const volScalarField& TCoolOrig, 
            const meshToMesh& neutroToFluid);  

            virtual void getFieldsLiquidFuel(
            const volVectorField& UOrig, 
            const volScalarField& porosityOrig, 
            const volScalarField& TfuelOrig, 
            const volScalarField& TcladOrig, 
            const volScalarField& rhoCoolOrig, 
            const volScalarField& TCoolOrig, 
            const meshToMesh& neutroToFluid,
            const rhoThermo& thermo,
            const compressible::turbulenceModel& turb);  

            //- Deform mesh
            virtual void deformMesh(const meshToMesh& TMToNeutro,const volVectorField& DispOrig);  

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
