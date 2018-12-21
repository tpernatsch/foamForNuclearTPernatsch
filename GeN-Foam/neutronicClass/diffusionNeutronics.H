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

        //Data
        IOdictionary nuclearData;
        IOdictionary nuclearDataRadialExp;
        IOdictionary nuclearDataAxialExp;
        IOdictionary nuclearDataFuelTemp;
        IOdictionary nuclearDataRhoCool;
        IOdictionary nuclearDataTCool;
        IOdictionary nuclearDataCladExp;
        label energyGroups;
        label precGroups;

        //Data fields
        PtrList<volScalarField> IV;
        PtrList<volScalarField> D;
        volScalarField Dalbedo;
        PtrList<volScalarField> nuSigmaEff;
        PtrList<volScalarField> sigmaPow;
        PtrList<volScalarField> sigmaDisapp;
        List<PtrList<volScalarField > > sigmaFromTo;
        PtrList<volScalarField> chiPrompt;
        PtrList<volScalarField> chiDelayed;
        PtrList<volScalarField> Beta;
        volScalarField BetaTot;
        PtrList<volScalarField> lambda;
        volScalarField fuelFraction;
        PtrList<volScalarField > discFactor;

        //Fields
        PtrList<volScalarField> flux;
        PtrList<volScalarField> fluxStar;
        PtrList<volScalarField> prec;
        volScalarField fluxStarAlbedo;
        volScalarField defaultFlux;
        volScalarField defaultPrec;

        //Derived fields
        volScalarField oneGroupFlux;
        volScalarField neutroSource;
        volScalarField delayedNeutroSource;
        volScalarField scatteringSourceExtra;
        volScalarField radExp;
        volScalarField axExp;

        //mappedFields
        volVectorField U;
        surfaceScalarField phi;
        volScalarField porosity;
        volScalarField Tfuel;
        volScalarField Tclad;
        volScalarField rhoCool;
        volScalarField TCool;
        volScalarField diffCoeffPrec;

        //XS
            label zoneNumber;

        //Disc factors
            // fuel fraction used to determine vol power IN the fuel starting from vol power from sigmaPow
            PtrList<scalar > fuelFractionList;
            // adjust or not DF
            PtrList<bool > dfAdjustList;
            //discontinuity factors
            PtrList<scalarField > discFactorList;
            //integral fluxes for adapting disc factors
            PtrList<scalarField > integralFluxList;

        //Feedbacks Coefficients

            bool fastNeutrons;
            bool adjustDiscFactors;
            bool useGivenDiscFactors;
            List<int> groupsWoDF;
            List<int> doNotParametrize;
            PtrList<entry> entries;
            //Constants
            PtrList<scalarField > IVList;
            PtrList<scalarField > chiPromptList;
            PtrList<scalarField > chiDelayedList;
            PtrList<scalarField > BetaList;
            PtrList<scalar > BetaTotList;
            PtrList<scalarField > lambdaList;
            //Nominal 
            PtrList<scalarField > DList;
            PtrList<scalarField > nuSigmaEffList;
            PtrList<scalarField > sigmaPowList;
            PtrList<scalarField > sigmaDisappList;
            PtrList<scalarSquareMatrix > sigmaFromToList; 
            // fuel temperature coefficient
            scalar TfuelRef;
            scalar TfuelPerturbed;
            PtrList<scalarField > fuelTempDList;
            PtrList<scalarField > fuelTempNuSigmaEffList;
            PtrList<scalarField > fuelTempSigmaPowList;
            PtrList<scalarField > fuelTempSigmaDisappList;
            PtrList<scalarSquareMatrix > fuelTempSigmaFromToList;
            // axial expansion coefficient
            scalar AxExp;
            PtrList<scalarField > axialExpDList;
            PtrList<scalarField > axialExpNuSigmaEffList;
            PtrList<scalarField > axialExpSigmaPowList;
            PtrList<scalarField > axialExpSigmaDisappList;
            PtrList<scalarSquareMatrix > axialExpSigmaFromToList;
            // radial expansion coefficient
            scalar RadExp;
            vector axialOrientation;
            PtrList<scalarField > radialExpDList;
            PtrList<scalarField > radialExpNuSigmaEffList;
            PtrList<scalarField > radialExpSigmaPowList;
            PtrList<scalarField > radialExpSigmaDisappList;
            PtrList<scalarSquareMatrix > radialExpSigmaFromToList;
            // coolant expansion coefficient
            scalar rhoCoolRef;
            scalar rhoCoolPerturbed;
            PtrList<scalarField > rhoCoolDList;
            PtrList<scalarField > rhoCoolNuSigmaEffList;
            PtrList<scalarField > rhoCoolSigmaPowList;
            PtrList<scalarField > rhoCoolSigmaDisappList;
            PtrList<scalarSquareMatrix > rhoCoolSigmaFromToList;
            // coolant temperature coefficient
            scalar TCoolRef;
            scalar TCoolPerturbed;
            PtrList<scalarField > TCoolDList;
            PtrList<scalarField > TCoolNuSigmaEffList;
            PtrList<scalarField > TCoolSigmaPowList;
            PtrList<scalarField > TCoolSigmaDisappList;
            PtrList<scalarSquareMatrix > TCoolSigmaFromToList;
            // cladding expansion coefficient
            scalar TcladRef;
            scalar TcladPerturbed;
            PtrList<scalarField > cladExpDList;
            PtrList<scalarField > cladExpNuSigmaEffList;
            PtrList<scalarField > cladExpSigmaPowList;
            PtrList<scalarField > cladExpSigmaDisappList;
            PtrList<scalarSquareMatrix > cladExpSigmaFromToList;

        //For control rods
        IOdictionary CRmove;
        PtrList<entry> CRentries;
        label CRNumber;
        PtrList<scalar > CRstart;
        PtrList<scalar > CRfinish;
        PtrList<scalar > CRspeed;
        PtrList<word >   CRFollowerName;
        PtrList<scalar > CRinitialPosition;
        PtrList<scalar > CRposition;
        PtrList<scalar > initialDistanceFromMeshCR;

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
                return flux;
            }

            //- Return precursors
            inline const PtrList<volScalarField> precursors() const
            {
                return prec;
            }

            //- Return one group flux
            inline const volScalarField& Flux() const
            {
                return oneGroupFlux;
            }

        // Calc

            //- Correct/update the properties
            virtual void correct(const label couplingIter, scalar& residual);  

            //- Take from another mesh the fields for XS parametrization
            virtual void getFields(
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
