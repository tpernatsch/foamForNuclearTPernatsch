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
#include "pimpleControl.H"
#include "saturationModel.H"
#include "typedefs.H" //- To access typedef of volScalarFieldTable

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;
class structureModel;

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
    const pimpleControl& pimple_;

    //- Ref to fluid1
    const fluid& fluid1_;

    //- Ref to fluid2
    const fluid& fluid2_;

    //- Ref to structure
    const structureModel& structure_;

    //- Ref to system pressure
    const volScalarField& p_;

    //- Ref to heat transfer coefficients table
    const volScalarFieldPtrTable& htcs_;

    //- Refs to volumetric mass transfer, interfacial temperature and area
    volScalarField& dmdt_;
    volScalarField& iT_;
    volScalarField& iA_;
   
    //- Ptr to saturation model
    autoPtr<saturationModel> saturation_;

    //- Minimum interfacial area
    const dimensionedScalar residualIA_;

    //- Cells in which the interfacial area is limited to a minimum so that
    //  boilng can start
    labelList residualIACells_;

    //- Total cumulative mass that should have been transferred between the 
    //  phases that was removed by the dmdt limiter 
    scalar dmLostToLimiter_;

public:

    TypeName("phaseChangeModel");

    //- Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        phaseChangeModel,
        phaseChangeModels,
        (
            const dictionary& dict,
            const pimpleControl& pimple, 
            const fluid& fluid1,
            const fluid& fluid2,
            const volScalarField& p,
            const volScalarFieldPtrTable& htcs,
            volScalarField& dmdt,
            volScalarField& iT,
            volScalarField& iA

        ),
        (dict, pimple, fluid1, fluid2, p, htcs, dmdt, iT, iA)
    );

    //- Constructors

        phaseChangeModel
        (
            const dictionary& dict,
            const pimpleControl& pimple, 
            const fluid& fluid1,
            const fluid& fluid2,
            const volScalarField& p,
            const volScalarFieldPtrTable& htcs,
            volScalarField& dmdt,
            volScalarField& iT,
            volScalarField& iA
        );

        //- Destructor
        virtual ~phaseChangeModel(){}

    //- Selectors

        static autoPtr<phaseChangeModel> New
        (
            const dictionary& dict,
            const pimpleControl& pimple, 
            const fluid& fluid1,
            const fluid& fluid2,
            const volScalarField& p,
            const volScalarFieldPtrTable& htcs,
            volScalarField& dmdt,
            volScalarField& iT,
            volScalarField& iA
        );

    //- Member Functions

        //- Set interfacial temperature to saturation temperature, under-relax
        //  if necessary
        void correctInterfacialT();

        //- Set interfacial area to a minimum so boiling can start (crude, good
        //  for now)
        void limitInterfacialArea();

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

        //- Calculate mass transfer term
        virtual void correctMassTransfer() = 0;

        //- Access

            const saturationModel& saturation() const
            {
                return saturation_();
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
