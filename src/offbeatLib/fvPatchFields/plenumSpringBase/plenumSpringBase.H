/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::plenumSpringBase

Description
    Base class with common functinality for plenumSpring patch fields.

SourceFiles
    plenumSpringBase.C

\*---------------------------------------------------------------------------*/

#ifndef plenumSpringBase_H
#define plenumSpringBase_H

#include "fvMesh.H"
#include "volFields.H"
#include "PrimitivePatchInterpolation.H"
#include "volPointInterpolation.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class plenumSpringBase Declaration
\*---------------------------------------------------------------------------*/

class plenumSpringBase
{
private:

    // Private data

        //- Patch
        const fvPatch& patch_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        plenumSpringBase(const plenumSpringBase&);

        //- Disallow default bitwise assignment
        void operator=(const plenumSpringBase&);

protected:

    // Protected data

        //- Spring modulus
        scalar springModulus_;   

        //- Initial spring loading (positive for compression, in m)
        scalar springPreCompression_;   

        //- Name list of top fuel surfaces.
        wordList fuelTopPatchNames_;

        //- IDs of top fuel surfaces.
        labelList fuelTopPatchIDs_;

        //- Name list of top cap inner surfaces. If the cap is not explicitly 
        // modeled, select the top cladding surface
        wordList topCapInnerPatchNames_;

        //- ID of patches connected to the top of the spring.
        labelList topCapInnerPatchIDs_;

        //- Name of strain field
        word displacementName_;  

public:

    //- Runtime type information
    TypeName("plenumSpringBase");


    // Constructors

        //- Construct from patch
        plenumSpringBase(const fvPatch&);

        //- Construct from dictionary
        plenumSpringBase(const fvPatch&, const dictionary&);

         //- Construct as copy, resetting patch
        plenumSpringBase(const fvPatch&, const plenumSpringBase&);


    //- Destructor
    virtual ~plenumSpringBase();


    // Member Functions


        // Access

            //- Return spring modulus
            virtual scalar springModulus() const
            {
                return springModulus_;
            };

            //- Return spring pre loading
            virtual scalar springPreCompression() const
            {
                return springPreCompression_;
            };


        // Operations

            //- Return total spring elongation
            virtual scalar springElongation() const;


        // IO

            //- Write
            virtual void write(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
