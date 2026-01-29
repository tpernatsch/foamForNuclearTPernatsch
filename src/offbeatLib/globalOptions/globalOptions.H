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
    Foam::globalOptions

Description
    This class reads global options from the main solver dictionary. It also
    calculates or read the angular fraction of the model.

    The options are assumed to be stored in a subDictionary named globalOptions.
    These options are made accessible with access functions.

SourceFiles
    globalOptions.C

\todo
    - globalOptions is necessary because some variables are needed by multiple 
    classe. However, the number of variables that require this treatement is
    limited for the moment. 
    Therefore, we should carefully consider if we can avoid having globalOptions 
    class at all.

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef globalOptions_H
#define globalOptions_H

#include "fvMesh.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class globalOptions Declaration
\*---------------------------------------------------------------------------*/

class globalOptions
:
    public regIOobject
{
    // Private data
    
        //- Reference to object fvMesh
        const fvMesh& mesh_;

        //- GlobalOptions dictionary
        const dictionary globalOptDict_;
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        globalOptions(const globalOptions&);

        //- Disallow default bitwise assignment
        void operator=(const globalOptions&);

protected:
    
    // Protected data  

        //- Reactor type (e.g. LWR)
        const word reactorType_;

        // Main axial direction
        vector pinDirection_;

        // Main radial direction
        vector radialDirection_;

        //- Angular fraction of modeled geometry 
        //- (e.g. 0.25 for a quarter of a pin, or 2/360 for a 2degree wedge )
        mutable scalar angularFraction_;

        void calcAngularFraction(const dictionary& dict) const;

public:

    // Constructors

        //- Construct from mesh
        globalOptions
        (
            const fvMesh& mesh,
            const dictionary& solverDict
        );
        

    //- Destructor
        ~globalOptions()
        {};


    //- Access functions
        const dictionary globalOptDict()
        {
            return globalOptDict_;
        }

        const word reactorType() const
        {
            return reactorType_;
        }

        const vector pinDirection() const
        {
            return pinDirection_;
        }

        const vector radialDirection() const
        {
            return radialDirection_;
        }

        const scalar& angularFraction() const
        {
            // Not set yet
            if(angularFraction_ == -1 or mesh_.topoChanging())
            {
                calcAngularFraction(globalOptDict_);
            }

            return angularFraction_;
        }


    // Member Functions
            

    // IO   
        
        //- Read data from time directory
        virtual bool readData(Istream& is)
        {
            return !is.bad();
        }
        
        //- Write data to time directory
        virtual bool writeData(Ostream& os) const 
        {
            return os.good();
        };    
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
