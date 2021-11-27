/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2018 OpenFOAM Foundation
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
    Foam::phaseChangeModel

Description

SourceFiles
    phaseChangeModel.C
    newphaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "IOdictionary.H"
#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "customPimpleControl.H"
#include "latentHeatModel.H"
#include "saturationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FFPair;

//- For the friend declaration (which is probably unnecessary)
namespace FSHeatTransferCoefficients
{
    class multiRegimeBoiling;
}

/*---------------------------------------------------------------------------*\
            Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
:
    public IOdictionary
{
protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to pimple dict in fvSolution to read under-relaxation factors
    const customPimpleControl& pimple_;

    //-
    FFPair& pair_;

    //- Ref to fluid1
    const fluid& fluid1_;

    //- Ref to fluid2
    const fluid& fluid2_;

    //- Ref to liquid
    const fluid& liquid_;

    //- Ref to vapour
    const fluid& vapour_;

    //- Ref to system pressure
    const volScalarField& p_;

    //- Refs to fluid temperatures
    const volScalarField& T1_;
    const volScalarField& T2_;

    //- Ref to interfacial heat transfer coefficients
    const volScalarField& htc1_;
    const volScalarField& htc2_;

    //- Refs to pair interfacialTemperature and area
    volScalarField& iT_;
    volScalarField& iA_;

    //- Latent heat
    volScalarField L_;

    //- Mass transfer (from fluid1 to fluid2) at the fluid-fluid interface
    volScalarField dmdtI_;

    //- Mass transfer (from fluid1 to fluid2) at the fluid1-wall (if > 0) or 
    //  fluid2-wall (if < 0) interface (e.g. sub-cooled boiling)
    volScalarField dmdtW_;

    //- Total mass transfer
    volScalarField dmdt_;

    //- Ptr to latentHeatModel
    autoPtr<latentHeatModel> latentHeatPtr_;
   
    //- Ptr to saturation model
    autoPtr<saturationModel> saturationPtr_;

    //-
    scalar residualIA_;

    //-
    DynamicList<label> residualIACells_;

    //- Total cumulative mass that should have been transferred between the 
    //  phases that was removed by the dmdt limiter 
    scalar dmLostToLimiter_;

    //- Linearized semi-implicit heat sources due to mass transfer (inclusive
    //  of intrinsic mass transfer enthalpy contribution)
    typedef HashTable<fvScalarMatrix*, word, word::hash> 
        fvScalarMatrixPtrTable;
    fvScalarMatrixPtrTable heSources_;

    //- Protected member functions

    //- Set interfacial area to a minimum so boiling can start (crude, good
    //  for now)
    void limitInterfacialArea();

    //- Set interfacial temperature to saturation temperature, under-relax
    //  if necessary
    void correctInterfacialTemperature();

    //- Limit mass transfer in two ways:
    //  -   make it so that mass cannot be removed from a fluid in a cell
    //      if there is no fluid to remove mass from in that cell (it is
    //  -   apply an adaptive mass transfer limiter to limit the impact of
    //      sporious continuity error contributions on the parameters that
    //      govern mass transfer (e.g. temperatures) in cells where the
    //      continuity error is excessive. The idea behind the limiter is
    //      described in the source file. Am I sure this is needed? Nope,
    //      but this PhD I'm doing is basically about testing whatever
    //      bullshit idea I can come up with to make things better. Jesus
    //      fucking Christ, I hate this
    void limitMassTransfer();

public:

    TypeName("phaseChangeModel");

    //- This was added by
    //  Carlo to make sure that multiRegimeBoiling can access dmdtW_ without
    //  issues. However, I believe this to be unnecessary as dmdtW_ is accessed
    //  via the objectRegistry lookup mechanism (in multiRegimeBoiling).
    //  If friend keywords were necessary, you would see them all over 
    //  OpenFOAM wherever objectRegistry lookups are involved, and in any
    //  class whose members (or the class itself) could be targets of a lookup.
    //  This is not the case
    friend class Foam::FSHeatTransferCoefficients::multiRegimeBoiling;

    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseChangeModel,
        phaseChangeModels,
        (
            FFPair& pair,
            const dictionary& dict
        ),
        (pair, dict)
    );

    //- Constructors

        phaseChangeModel
        (
            FFPair& pair,
            const dictionary& dict
        );

        //- Destructor
        virtual ~phaseChangeModel(){}

    //- Selectors

        static autoPtr<phaseChangeModel> New
        (
            FFPair& pair,
            const dictionary& dict
        );

    //- Member Functions

        //- Calculate mass transfer term
        virtual void correct();

        //- Calculate interfacial mass transfer term
        virtual void correctInterfacialDmdt() = 0;

        //- Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const FFPair& pair() const
            {
                return pair_;
            }

            const fluid& liquid() const
            {
                return liquid_;
            }

            const fluid& vapour() const
            {
                return vapour_;
            }

            const volScalarField& L() const
            {
                return L_;
            }

            const volScalarField& dmdt() const
            {
                return dmdt_;
            }

            const fvScalarMatrix& heSource(const word& heName)
            {
                return *heSources_[heName];
            }

            /*
            const volScalarField& L() const
            {
                return latentHeat_->L();
            }

            const saturationModel& saturation() const
            {
                return saturation_();
            }*/
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
