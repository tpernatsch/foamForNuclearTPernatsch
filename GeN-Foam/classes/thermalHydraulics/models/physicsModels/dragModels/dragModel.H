/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::dragModel

Description
    Class to handle interfacial drag between either the two fluid phases
    or a fluid phase and the structure. It is run-time selectable depending
    on whether a fluid-fluid (FF in short) or a fluid-structure (FS in short)
    pair is provided. For this reason, the class has two constructors
    and, consequently, two selector implementations in newDragModel.C. However,
    as the return type of Kd is always a tensorField, only one base class was
    deemed necessary, in order to keep it simple.
    Fluid-fluid models are to be added to the FF selection table, while 
    fluid-structure models to the FS selection table.

    By Stefan Radman (stefan.radman@epfl.ch)

SourceFiles
    dragModel.C
    newDragModel.C

\*---------------------------------------------------------------------------*/

#ifndef dragModel_H
#define dragModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluid.H"
#include "structureModel.H"
#include "FFPair.H"
#include "FSPair.H"
#include "zeroGradientFvPatchFields.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class dragModel Declaration
\*---------------------------------------------------------------------------*/

class dragModel
:
    public IOdictionary
{
public :

    //- Only relevant for fluid-structure drag
    enum transverseDragModel
    {
        isotropic,
        same,
        Blasius, 
        GunterShaw
    };
    static const Enum<transverseDragModel> transverseDragModelNames_;

protected:

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Name of phase pair
    const word pairName_;

    //- Ptr to FFPair if this dragModel is between fluid and fluid
    const FFPair* const FFPair_;

    //- Ptr to FSPair if this dragModel is between fluid and structure
    const FSPair* const FSPair_;

    //- List of cells in which the dragModel exists (only for FS  and shared
    //  dragModels)
    labelList cellList_;

    //- Field of 1s where in which the dragModel exists (only for FS 
    //  and shared dragModels)
    scalarField cellField_;

    //- The following is only for fluid-structure dragModels
        
        //-
        const transverseDragModel transverseDragModel_;

        //- Minuscule optimizations...
        bool isotropic_;

        //- These three labels tell me the orientation of the two transverse
        //  axes (dta0,1) and the principal axis (dpa_). The values range from 0 to
        //  2 and are 0 for localX, 1 for localY, 2 for localZ
        label dta0_;
        label dta1_;
        label dpa_;

        //- i.e. 0.5*alpha*rho*magU, not a volScalarField because the dragModel
        //  only affects certain cellZone(s). Alpha is the phase fraction of
        //  the (only) fluid in fluid-structure drag, rho its density, magU its
        //  velocity magnitude
        mutable scalarField* halfAlphaRhoMagU_;

    //- Protected member functions

        //- Correct drag components perpendicular to main flow direction. Only 
        //  used for fluid-structure drag (fluid-fluid drag in inherently 
        //  isotropic)
        void correctTransverseKd(volTensorField& Kd) const;

public:

    TypeName("dragModel");

    //- Static data members
    static const dimensionSet dimK;

    //- Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragModel,
            FFDragModels,
            (
                const objectRegistry& objReg,
                const dictionary& dict,
                const FFPair& FFPair
            ),
            (objReg, dict, FFPair)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            dragModel,
            FSDragModels,
            (
                const objectRegistry& objReg,
                const dictionary& dict,
                const FSPair& FSPair
            ),
            (objReg, dict, FSPair)
        );

    //- Constructors

        //- Fluid-fluid drag
        dragModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair
        );

        //- Fluid-structure drag
        dragModel
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair
        );

        //- Destructor
        virtual ~dragModel(){}

    //- Selectors

        static autoPtr<dragModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FFPair& FFPair
        );

        static autoPtr<dragModel> New
        (
            const objectRegistry& objReg,
            const dictionary& dict,
            const FSPair& FSPair
        );

        //- Make in table given a key, dict, source and destination object
        //  registries. Source registers are those where the FF/FS pairs will
        //  be looked for. Destination ones are those in which the model 
        //  IOobject will be registered
        static void makeInTable
        (
            const objectRegistry& srcObjReg,
            const objectRegistry& dstObjReg,
            const word key,
            const dictionary& dict,
            HashTable<autoPtr<dragModel>, word, word::hash>& table
        );

    //- Member Functions

        //- Drag coefficient used in momentum equations:
        //    ddt(alpha1*rho1*U1) + ... = ... Kd&(U2-U1)
        //    ddt(alpha2*rho2*U2) + ... = ... Kd&(U1-U2)
        //  for fluid-fluid drag or 
        //    ddt(alpha_i*rho_i*U_i) + ... = ... -Kd&U_i
        //  for fluid-structure drag.
        //  Kd is a tensor for fluid-structure drag and a
        //  scalar for fluid-fluid drag, but for programming
        //  convenience, fluid-fluid Kd is also represented
        //  as a diagonal tensor with all equal entries.
        virtual void correctKd(volTensorField& Kd) const = 0;

        //- Access

            const word pairName() const
            {
                return pairName_;
            }

            const labelList& cellList() const
            {
                return cellList_;
            }

            const scalarField& cellField() const
            {
                return cellField_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
