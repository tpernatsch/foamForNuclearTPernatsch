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
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::heatSource

Description
    Parent class for heat source model. 

    \htmlinclude heatSource_heatSource.html
    
SourceFiles
    heatSource.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)\n

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef heatSource_H
#define heatSource_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "materials.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heatSource Declaration
\*---------------------------------------------------------------------------*/

class heatSource
{
    // Private data

    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        heatSource(const heatSource&);

        //- Disallow default bitwise assignment
        void operator=(const heatSource&);
        
 protected:   

    // Protected data
        
        //- Reference to mesh
        const fvMesh& mesh_;

        //- Reference to the materials class
        const materials& mat_;

        //- Reference to the heatSourceOptions dictionary
        const dictionary heatSourceOptDict_;
        
        //- Heat source field [W/m3]
        volScalarField Q_;

        //- Value of current time
        scalar currentTime_;      

        //- Maximum relative power increase for adjustable time step
        scalar maxRelativePowerIncrease_;     

        //- Maximum relative power decrease for adjustable time step
        scalar maxRelativePowerDecrease_;   

        //- Maximum absolute power increase for adjustable time step
        scalar maxAbsolutePowerIncrease_;

    
    // Protected Member Functions

public:

    //- Runtime type information
    TypeName("fromLatestTime");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            heatSource,
            dictionary,
            ( 
                const fvMesh& mesh,
                const materials& materials,
                const dictionary& heatSourceOptDict
            ),
            (mesh, materials, heatSourceOptDict)
        );


    // Constructors

        //- Construct from mesh, materials and dict
        heatSource
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceOptDict
        );


    // Selectors

        //- Select from dictionary
        static autoPtr<heatSource> New
        (        
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& solverDict
        );


    //- Destructor
    virtual ~heatSource();


    // Access functions
      
        //- Return non-cost reference to heat source
        virtual volScalarField& Q()
        {
            return Q_;
        }  
          
        //- Return cost reference to heat source      
        virtual const volScalarField& Q() const
        {
            return Q_;
        }  
      
        //- Return cost reference to mesh
        virtual const fvMesh& mesh() const
        {
            return mesh_;
        } 

        
    // Member Functions


        //- Update power distribution
        virtual void correct(); 

        //- Return the average power density 
        virtual scalar averagePowerDensity() const
        {
            return GREAT;
        }; 
            
        //- Return the next average power density (assuming the deltaT stays the
        // same)
        virtual scalar predictNextPowerDensity() const
        {
            return GREAT;
        };          

        //- Provide next deltaT due to absolute change in heat source (not relative)
        virtual scalar nextDeltaTAbsoluteChange(scalar deltaT) const;
          
        //- Return the final time marker after which the power distribution
        // is no longer defined.
        virtual scalar lastTimeMarker() const
        {
            return GREAT;
        };        

        //- For time-dependent heat sources, return the next time marker,
        // i.e. the next time point at which the rate of power change
        // will be discontinuous or non-smooth.
        virtual scalar nextTimeMarker() const
        {
            return GREAT;
        }
              
        //- Return next allowable time step based on rate of change of power
        virtual scalar nextDeltaT();  

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
