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
    Foam::thermoMechanicsSolver

Description
Class for XS parametrization


\*---------------------------------------------------------------------------*/

#ifndef XS_H
#define XS_H

#include "volFields.H"
#include "rhoThermo.H"
#include "fvMatrices.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "coordinateSystem.H"
#include "fvCFD.H"
#include "fvMatrixExt.H"
#include "turbulentFluidThermoModel.H"
#include "volPointInterpolation.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                         Class thermoMechanicsSolver Declaration
\*---------------------------------------------------------------------------*/

class XS
{
    // Private data

    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        XS(const XS&);

        //- Disallow default bitwise assignment
        void operator=(const XS&);
        
protected:
    
    // Protected data

        //- Reference to mesh
            const fvMesh& mesh_;

        //Data
            IOdictionary nuclearData;
            IOdictionary nuclearDataRadialExp;
            IOdictionary nuclearDataAxialExp;
            IOdictionary nuclearDataFuelTemp;
            IOdictionary nuclearDataRhoCool;
            IOdictionary nuclearDataTCool;
            IOdictionary nuclearDataCladExp;

            label energyGroups_;
            label precGroups_;
        
        //Data fields
            PtrList<volScalarField> IV_;
            PtrList<volScalarField> D_;
            PtrList<volScalarField> nuSigmaEff_;
            PtrList<volScalarField> sigmaPow_;
            PtrList<volScalarField> sigmaDisapp_;
            List<PtrList<volScalarField > > sigmaFromTo_;
            PtrList<volScalarField> chiPrompt_;
            PtrList<volScalarField> chiDelayed_;
            PtrList<volScalarField> Beta_;
            volScalarField BetaTot_;
            PtrList<volScalarField> lambda_;
            volScalarField fuelFraction_;
            PtrList<volScalarField > discFactor_;
            scalar ScNo_;
        
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
            bool adjustDiscFactors_;
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

public:

    //- Runtime type information
    TypeName("XS");
    
    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        XS,
        dictionary,
        (const fvMesh& mesh),
        (mesh)
    );    
    // Access

        inline const label energyGroups() const
        {
            return energyGroups_;
        }
        inline const label precGroups() const
        {
            return precGroups_;
        }
        inline const PtrList<volScalarField> IV() const
        {
            return IV_;
        }
        inline const PtrList<volScalarField> D() const
        {
            return D_;
        }
        inline const PtrList<volScalarField> nuSigmaEff() const
        {
            return nuSigmaEff_;
        }
        inline const PtrList<volScalarField> sigmaPow() const
        {
            return sigmaPow_;
        }
        inline const PtrList<volScalarField> sigmaDisapp() const
        {
            return sigmaDisapp_;
        }
        inline const List<PtrList<volScalarField > > sigmaFromTo() const
        {
            return sigmaFromTo_;
        }
        inline const PtrList<volScalarField> chiPrompt() const
        {
            return chiPrompt_;
        }
        inline const PtrList<volScalarField> chiDelayed() const
        {
            return chiDelayed_;
        }
        inline const PtrList<volScalarField> Beta() const
        {
            return Beta_;
        }
        inline const volScalarField BetaTot() const
        {
            return BetaTot_;
        }
        inline const PtrList<volScalarField> lambda() const
        {
            return lambda_;
        }
        inline const volScalarField fuelFraction() const
        {
            return fuelFraction_;
        }
        inline const PtrList<volScalarField> discFactor() const
        {
            return discFactor_;
        }
        inline const scalar ScNo() const
        {
            return ScNo_;
        }

    // Constructors

        //- Construct from mesh and thermo-mechanical properties
        XS
        (
            const fvMesh& mesh
        );
        
        //- Return a pointer to a new XS
        static autoPtr<XS> New
        (
            const fvMesh& mesh
        );           

    //- Destructor
        virtual ~XS();

    // Member Functions
    
        //- Update XS      

        void correct
        (
            const volScalarField& Tfuel, 
            const volScalarField& Tclad, 
            const volScalarField& rhoCool, 
            const volScalarField& TCool,
            const volVectorField& Disp             
        );

        void init();

        void adjustDiscFactors(const PtrList<volScalarField>& fluxStar);


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //