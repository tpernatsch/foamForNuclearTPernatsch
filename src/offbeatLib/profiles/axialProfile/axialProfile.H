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
    Foam::axialProfile

Description
    Base class for a general profile in the axial direction. The base
    class assumes a uniform "flat" profile.
    The class is provided with a mesh and addressing array, and then deals with
    the mapping of the profile onto an addressed list, which can be accessed
    through the profile function.

Usage
    \verbatim
    {
        type    "flat";
    }
    \endverbatim
    
SourceFiles
    axialProfile.C

\mainauthor
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef axialProfile_H
#define axialProfile_H

#include "primitiveFields.H"
#include "objectRegistry.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class axialProfile Declaration
\*---------------------------------------------------------------------------*/

class axialProfile
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        axialProfile(const axialProfile&);

        //- Disallow default bitwise assignment
        void operator=(const axialProfile&);

protected:
    
    // Protected data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Addressing list
        const labelList& addr_;

        //- Profile
        scalarField profile_;

        //- Scaling factor profile
        scalarField scalingFactorProfile_;
        
        //- Axial direction
        vector axialDirection_;
        
        //- Minimum axial location
        scalar zMin_;
        
        //- Maximum axial location
        scalar zMax_;

public:

    //- Runtime type information
    TypeName("flat");

    // Declare run-time constructor selection table

    declareRunTimeSelectionTable
    (
        autoPtr,
        axialProfile,
        dictionary,
        (
            const fvMesh& mesh, 
            const dictionary& dict, 
            const labelList& addr,
            const vector& axialDirection,
            scalar zMin,
            scalar zMax 
        ),
        (mesh, dict, addr, axialDirection, zMin, zMax)
    );

    // Constructors

        //- Construct from dictionary
        axialProfile
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            scalar zMin,
            scalar zMax 
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<axialProfile> New
        (
            const fvMesh& mesh,    
            const dictionary& dict,
            const labelList& addr,
            const vector& axialDirection,
            scalar zMin,
            scalar zMax 
        );

    //- Destructor
    virtual ~axialProfile();


    // Member Functions
    
        //- Return a const reference to the mesh
        inline const fvMesh& mesh() const
        {
            return mesh_;
        };  

        //- Update the profile
        virtual void correct()
        {};

        //- Return a reference to the profile
        inline const scalarField& profile() const
        {
            return profile_;
        };

        //- Return a reference to the profile
        inline const scalarField& scalingFactorProfile() const
        {
            return scalingFactorProfile_;
        };
        
        //- Return the addressed axial locations as a field normalised
        //- between 0 and 1
        tmp<scalarField> axialLocations() const;
        
        //- For time-dependent heat sources, return the next time marker,
        //- i.e. the next time point at which the rate of change
        //- will be discontinuous or non-smooth.
        virtual scalar nextTimeMarker() const
        {
            return GREAT;
        }

        //- Return the final time marker after which the field distribution
        //- is no longer defined.
        virtual scalar lastTimeMarker() const
        {
            return GREAT;
        }

        //- Check that the axial profile is normalized to 1 with respect 
        //- to the volume
        virtual void checkNormalization() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
