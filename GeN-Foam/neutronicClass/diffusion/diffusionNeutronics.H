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
    Foam::neutronics::diffusionNeutronics

Description
    subclass of neutronics. Created to evaluate pressure drops and
    heat transfer in a porous medium using available correlations.
    Different correlations can be used for three principal directions
    of a local coordinate system defined by the user. The possibility
    also exists to use different correlations for laminar and turbulent
    flow.

SourceFiles
    diffusionNeutronics.C

Reference publications
    Carlo Fiorina, Ivor Clifford, Manuele Aufiero, Konstantin Mikityuk,
    "GeN-Foam: a novel OpenFOAM® based multi-physics solver for 2D/3D transient
    analysis of nuclear reactors", Nuclear Engineering and Design, submitted

    Carlo Fiorina, Konstantin Mikityuk, " Application of the new GeN-Foam multi-physics
    solver to the European Sodium Fast Reactor and verification against available codes",
    Proceedings of ICAPP 2015, May 03-06, 2015 - Nice (France), Paper 15226


Author
    Carlo Fiorina <carlo.fiorina@outlook.com; carlo.fiorina@epfl.ch;>

\*---------------------------------------------------------------------------*/

#ifndef diffusionNeutronics_H
#define diffusionNeutronics_H

#include "neutronics.H"
#include "XS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class diffusionNeutronics Declaration
\*---------------------------------------------------------------------------*/

class diffusionNeutronics
:
    public neutronics//class inheritance
{
protected:

    // Protected data

        // Cross sections
        XS xs_;

        //Fields
        PtrList<volScalarField> flux_;
        PtrList<volScalarField> fluxStar_;
        PtrList<volScalarField> prec_;
        volScalarField fluxStarAlbedo_;
        volScalarField defaultFlux_;
        volScalarField defaultPrec_;
        volScalarField Dalbedo_;

        //Derived fields
        volScalarField oneGroupFlux_;
        volScalarField neutroSource_;
        volScalarField delayedNeutroSource_;
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


private:

    // Private data

    // Private Member Functions

        //- Disallow default bitwise copy construct
        diffusionNeutronics(const diffusionNeutronics&);

        //- Disallow default bitwise assignment
        void operator=(const diffusionNeutronics&);


public:

    //- Runtime type information
        TypeName("diffusionNeutronics");


    // Constructors

        //- Construct from mesh, read data from IOdictionaries
        diffusionNeutronics(fvMesh& mesh);

    //- Destructor
        virtual ~diffusionNeutronics();


    // Member Functions

        // Access

            //- Return fluxes
            inline const PtrList<volScalarField> fluxes() const
            {
                return flux_;
            }

            //- Return precursors
            inline const PtrList<volScalarField> precursors() const
            {
                return prec_;
            }

            //- Return one group flux
            inline const volScalarField& Flux() const
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
