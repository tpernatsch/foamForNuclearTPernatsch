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
    public scalarField,
    public objectRegistry
{
private :
    
    //- Ref to mesh
    const fvMesh& mesh_;

    //- Dictionary of regime
    const dictionary dict_;

    //- Regime name
    const word name_;

    //- True if this regime exists in at least one cell in the domain
    bool isCurrentlyPresent_;

    //- True if this regime is interpolated between two other regimes
    const bool isInterpolated_;

    //- Ref to regime1 if this regime is interpolated (otherwise set to *this)
    const regime& regime1_;

    //- Ref to regime2 if this regime is interpolated (otherwise set to *this)
    const regime& regime2_;

    //- Interpolation coefficients
        //- Only meaningful for interpolated regimes. It is a field btw. 0 and 1. 
        //  It weighs the contribution of regime1 to the total regime.
    scalarField coeffs1_;

        //- Only meaningful for interpolated regimes. It is a field btw. 0 and 1. 
        //  It weighs the contribution of regime2 to the total regime.
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
            const regime& regime1,
            const regime& regime2
        );

        //- Destructor
        virtual ~regime(){}

    //- Members

        //- Adjust the isCurrentlyPresent_ flag according to the new spatial
        //  regime distribution
        void correct();

        //- Update the drag coefficients in a table according to the
        //  spatial regime distribution
        void correctDragTable(volTensorFieldTable& table) const;

        //- Update the heat transfer coefficients in a table according to
        //  the spatial regime distribution
        void correctHeatTransferTable(volScalarFieldTable& htcTable) const;

        //- Update the:
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
        void correctFluidGeometry
        (
            volScalarField& iA,
            volScalarField& frac1,
            volScalarField& frac2,
            fluid& fluid1,
            fluid& fluid2
        ) const;

        //- Access

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
                return isCurrentlyPresent_;
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

        //- Operators

            //- Make the = and [] operators of regime behave as those of
            //  its scalarField base class
            
            template<class T>
            inline void operator=(T& rhs)
            {
                this->scalarField::operator=(rhs);
            }

            template<class T>
            inline void operator=(tmp<T> rhs)
            {
                this->scalarField::operator=(rhs);
            }
            
            inline scalar& operator[](int celli)
            {
                return this->scalarField::operator[](celli);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
