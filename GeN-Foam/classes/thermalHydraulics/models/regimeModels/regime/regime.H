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
    Foam::regime

Description
    Class representing a regime. The class itself is a scalarField (not a 
    volScalarField as I don't need boundary conditions) representing the
    presence/absence of the regime for each mesh cell. Class members consist
    of :
    -   regime name;
    -   tables of regimePhysicsModels (indexed by for each phase pair);
    -   a flag to specify wether the regime is interpolated or not;
    -   a flag to specify wether the regime is currently present in the domain;
    -   references to the regimes between which this regime has to be
        (linearly) interpolated.
    If the regime is interpolated, the value of the scalarField is the
    interpolation coefficient between the two regimes.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    regimeMapModel.C

\*---------------------------------------------------------------------------*/

#ifndef regime_H
#define regime_H

#include "volFields.H"
#include "dictionary.H"
#include "HashPtrTable.H"
#include "typedefs.H"
#include "dragModel.H"
#include "twoPhaseDragMultiplierModel.H"
#include "heatTransferModel.H"
#include "fluidGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class FFPair;
class fluid;
class structureModel;

/*---------------------------------------------------------------------------*\
                           Class regimeMapModel Declaration
\*---------------------------------------------------------------------------*/

class regime
:
    public objectRegistry
{
private :
    
    //- Ref to mesh
    const fvMesh& mesh_;

    //- Dictionary of regime
    const dictionary dict_;

    //- Regime name
    const word name_;

    bool requiresModelCorrection_;

    //- True if this regime is interpolated between two other regimes
    const bool isInterpolated_;

    //- Const ptr to regime1 if this regime is interpolated
    regime* const regime1_;

    //- Const ptr to regime2 if this regime is interpolated
    regime* const regime2_;

    //- Field of 1 or 0s that mark where the field is present
    scalarField cellField_;

    //- List of cells where the field is present
    labelList cellList_;

    //- Interpolation coefficients

    //- Only meaningful for interpolated regimes. It is a field btw. 0 and 1. 
    //  It weighs the contribution of regime1 to the total regime
    scalarField coeffs1_;

    //- Only meaningful for interpolated regimes. It is a field btw. 0 and 1. 
    //  It weighs the contribution of regime2 to the total regime
    scalarField coeffs2_;

    //- Models

    //- Drag model table
    typedef HashTable
    <
        autoPtr<dragModel>,
        word,
        word::hash
    > dragModelTable;
    dragModelTable dragModels_;

    //- Two phase drag multiplier model
    autoPtr<twoPhaseDragMultiplierModel> twoPhaseDragMultiplier_;

    //- Heat transfer model table
    typedef HashTable
    <
        autoPtr<heatTransferModel>,
        word,
        word::hash
    > heatTransferModelTable;
    heatTransferModelTable heatTransferModels_;

    //- Fluid geometry
    autoPtr<fluidGeometry> fluidGeometry_;

    //- Local drag coeff table
    volTensorFieldPtrTable Kds_;

    //- Local heat transfer coeff table
    volScalarFieldPtrTable htcs_;

public:

    //- Constructors

        //- For non-interpolated regimes (there is a dict)
        regime
        (
            const fvMesh& mesh,
            const word name,
            const dictionary& dict
        );

        //- For interpolated regimes (needless to say, regime1 and regime2 must
        //  have been already instantiated). The name of the interpolated
        //  regime will be regime1.name()+"."+regime2.name()
        regime
        (
            const fvMesh& mesh,
            regime& regime1,
            regime& regime2
        );

        //- Destructor
        virtual ~regime(){}

    //- Members

        //- Add keys found in dragTable, htcTable to Kds_ and htcs_
        //  respectively
        void initLocalTables
        (
            volTensorFieldPtrTable& KdTable, 
            volScalarFieldPtrTable& htcTable
        );

        //- Update local regime dragModels
        void correctDragModels();

        //- Update local regime heatTransferModel
        void correctHeatTransferModels();

        //- Update local regime fluidGeometryModels
        void correctFluidGeometryModels();

        //- Update the global drag coefficients in a table according to the
        //  spatial regime distribution
        void correctDragTable(volTensorFieldPtrTable& KdTable) const;

        //- Update the global heat transfer coefficients in a table according 
        //  to the spatial regime distribution
        void correctHeatTransferTable(volScalarFieldPtrTable& htcTable) const;

        //- Update the provided global:
        //      - fluid-fluid interfacial area;
        //      - fraction of structure interfacial area that contacts fluid1;
        //      - fraction of structure interfacial area that contacts fluid2;
        //      - characteristic dimension of fluid1;
        //      - characteristic dimension of fluid2;
        //      - marker field that specifies where fluid1 is continuous and
        //        where it is dispersed, with values between 0 and 1 in 
        //        interpolation regions;
        //      - marker field that specifies where fluid2 is continuous and
        //        where it is dispersed, with values between 0 and 1 in
        //        interpolation regions;
        // according to the spatial regime distribution
        void correctFluidGeometryFields
        (
            volScalarField& iA,
            volScalarField& frac1,
            volScalarField& frac2,
            fluid& fluid1,
            fluid& fluid2
        ) const;

        //- Access (const and non-const)

            const word name() const
            {
                return name_;
            } 

            bool isInterpolated() const
            {
                return isInterpolated_;
            }

            bool isCurrentlyPresent() const
            {
                return (cellList_.size() != 0);
            }

            bool requiresModelCorrection() const
            {
                return requiresModelCorrection_;
            }

            bool& requiresModelCorrection()
            {
                return requiresModelCorrection_;
            }

            const regime& regime1() const
            {
                return *regime1_;
            }

            const regime& regime2() const
            {
                return *regime2_;
            }

            regime& regime1()
            {
                return *regime1_;
            }

            regime& regime2()
            {
                return *regime2_;
            }

            const scalarField& cellField() const
            {
                return cellField_;
            }

            const labelList& cellList() const
            {
                return cellList_;
            }

            scalarField& cellField()
            {
                return cellField_;
            }

            labelList& cellList()
            {
                return cellList_;
            }

            const dragModelTable& dragModels() const
            {
                return dragModels_;
            }

            const heatTransferModelTable& heatTransferModels() const
            {
                return heatTransferModels_;
            }

            const fluidGeometry& fluidGeometryRef() const
            {
                return fluidGeometry_();
            }

            scalarField& coeffs1()
            {
                return coeffs1_;
            }

            scalarField& coeffs2()
            {
                return coeffs2_;
            }

            //- Local drag table
            const volTensorFieldPtrTable& Kds() const
            {
                return Kds_;
            }

            //- Local heat transfer table
            const volScalarFieldPtrTable& htcs() const
            {
                return htcs_;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
