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
    Foam::powerModel

Description
    Run-time selectable class to handle the thermal description (i.e. energy 
    source/sinks) of a structure.
    It returns a heat source term to the main solver be used in the energy 
    equations.

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

SourceFiles
    powerModel.C

\*---------------------------------------------------------------------------*/

#ifndef powerModel_H
#define powerModel_H

#include "dictionary.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "zeroGradientFvPatchFields.H"
#include "typedefs.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Forward declarations
class fluid;
class structureModel;

/*---------------------------------------------------------------------------*\
                           Class powerModel Declaration
\*---------------------------------------------------------------------------*/

class powerModel
:
    public IOdictionary
{
protected:

    //- Ref to parent structure
    const structureModel& structure_;

    //- Reference to mesh
    const fvMesh& mesh_;

    //- List of cells in which this powerModel exists
    labelList cellList_;

    //- Marker field for this power model, = 1 in all the cells that are
    //  listed in cellsList, and 0 elsewhere
    scalarField cellField_;

    //- Interfacial area of the powerModel
    volScalarField iA_;

    //- Volume fraction of the powerModel (can in principle be different from
    //  the structure alpha as not the entirety of a structure region might
    //  be producing power)
    volScalarField alpha_;

public:

    //- Runtime type information
    TypeName("powerModel");

    //- Runtime selection table declaration
    declareRunTimeSelectionTable
    (
        autoPtr,
        powerModel,
        powerModels,
        (
            const structureModel& structure,
            const dictionary& dicts
        ),
        (structure, dicts)
    );

    // Constructors

        //- Construct over specific cells
        powerModel
        (
            const structureModel& structure,
            const dictionary& dicts
        );

        //- Destructor
        virtual ~powerModel(){}

    // Selectors

        static autoPtr<powerModel> New
        (
            const structureModel& structure,
            const dictionary& dicts
        );

    // Members

        /*
        volScalarField& initFieldInTable
        (
            word name, 
            dimensionedScalar value,
            bool readIfPresent,
            bool autoWrite
        );
        */

        //- In two-phase solvers :
        //  HT = H1*frac1*T1 + H2*frac2*T2
        //  H =  H1*frac1 + H2*frac2
        //- In a mono-phase solver :
        //  HT = H*T
        //  H = H
        virtual void correct
        ( 
            const volScalarField& HT, 
            const volScalarField& H
        ) = 0;

        //- Set provided field with structure surface temperature
        virtual void correctT(volScalarField& T) const = 0;

        //- Turn off power
        virtual void powerOff(){}

        //- Remove abstractness inherited from regIOobject
        virtual bool writeData(Ostream& os) const 
        {
            return os.good();
        }

        // Access

            const volScalarField& iA() const
            {
                return iA_;
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
