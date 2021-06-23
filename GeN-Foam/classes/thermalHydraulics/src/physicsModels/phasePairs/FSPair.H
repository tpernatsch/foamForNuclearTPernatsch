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

Description
    Class to handle dimensionless numbers that characterise the system, e.g.
    Reynolds, etc. While some can (rightfully) think that some of these numbers
    should be directly members of other classes (e.g., Reynolds of fluids
    should be in the fluid class), in general these numbers depend on a
    multitude of parameters. You cannot define a Reynolds for a fluid by
    itself, it needs to be geometrically constratined in a certain domain 
    (so that a hydraulic diameter can be defined), so it is not properly
    "just" of the fluid. While this logic can of course be debate, I wanted
    to keep it simple and not start dispersing dimensionless numbers randomly
    around the code. 

Author
    Stefan Radman (stefan.radman@epfl.ch / stefanradman92@gmail.com)

Class
    Foam::FSPair

SourceFiles
    FSPair.C

\*---------------------------------------------------------------------------*/

#ifndef FSPair_H
#define FSPair_H

#include "fluid.H"
#include "structure.H"
#include "FSDragFactor.H"
#include "FSHeatTransferCoefficientModel.H"
#include "contactPartitionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class FSPair Declaration
\*---------------------------------------------------------------------------*/

class FSPair
    : public regIOobject
{
protected :

    //- Ref to mesh
    const fvMesh& mesh_;

    //- Ref to phaseProperties dict
    const dictionary& dict_;

    //-
    const dictionary& pimple_;

    //- Ref to fluid
    fluid& fluid_;

    //- Ref to structure
    structure& structure_; 

    //- Minimum allowable Reynolds
    const dimensionedScalar minRe_;

    //- Reynolds of the pair
    volScalarField Re_;

    //- local components of the Reynolds number in the local frame of
    //  reference directions
    volVectorField lRe_;

    //- Table of drag factors (per region name)
    typedef HashTable
    <
        autoPtr<FSDragFactor>,
        word,
        word::hash
    > KdTable;
    KdTable KdPtrs_;

    //- Total drag factor
    volTensorField Kd_;

    //-
    typedef HashTable
    <
        autoPtr<FSHeatTransferCoefficientModel>,
        word,
        word::hash
    > htcTable;
    htcTable htcPtrs_;

    //- Total heat transfer coefficient
    volScalarField htc_;

    //- Fluid-structure contact fraction
    autoPtr<contactPartitionModel> contactPartition_;
    autoPtr<volScalarField> fPtr_;

public:

    TypeName("FSPair");

    //- Constructors

        FSPair
        (
            fluid& fluid,
            structure& structure,
            const dictionary& dict
        );

        //- Destructor
        virtual ~FSPair(){}

    //- Member Functions

        //- To remove abstraction from regIOobject
        virtual bool writeData(Ostream& os) const
        {
            return os.good();
        }

        //-
        void correct
        (
            const bool& correctFluidDynamics, 
            const bool& correctEnergy
        );

        //- Access

            const fvMesh& mesh() const
            {
                return mesh_;
            }

            const word& fluidName() const
            {
                return fluid_.name();
            }

            fluid& fluidRef()
            {
                return fluid_;
            }

            const fluid& fluidRef() const
            {
                return fluid_;
            }

            structure& structureRef()
            {
                return structure_;
            }

            const structure& structureRef() const
            {
                return structure_;
            }

            word reversedName() const
            {
                return word("structure."+fluid_.name());
            }

            bool foundUnorderedPairSubDict
            (
                const dictionary& dict
            ) const
            {
                return 
                    (
                        dict.found(this->name())
                    or  dict.found(this->reversedName())
                    );
            }

            const dictionary& getUnorderedPairSubDict
            (
                const dictionary& dict
            ) const
            {
                return 
                    (dict.found(this->name())) ?
                    dict.subDict(this->name()) :
                    dict.subDict(this->reversedName());
            }

            const volScalarField& Re() const
            {
                return Re_;
            }

            const volVectorField& lRe() const
            {
                return lRe_;
            }    

            const volTensorField& Kd() const 
            {
                return Kd_;
            }

            volTensorField& Kd()
            {
                return Kd_;
            }

            const volScalarField& htc() const 
            {
                return htc_;
            }

            volScalarField& htc() 
            {
                return htc_;
            }

            autoPtr<volScalarField>& fPtr()
            {
                return fPtr_;
            }

            const autoPtr<volScalarField>& fPtr() const
            {
                return fPtr_;
            }

            const volScalarField& f() const
            {
                return fPtr_();
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
