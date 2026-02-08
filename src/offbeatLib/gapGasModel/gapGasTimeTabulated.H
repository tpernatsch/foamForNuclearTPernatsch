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
    Foam::gapGasTimeTabulated

Description
    Model for time-dependent gap gas composition, and pressure
    when the data is provided from external source (for instance, other fuel performance 
    codes). The time-dependent lists of mass fractions and gas pressures are provided 
    within the gapGasOptions subdictionary within the solverDict file.
    
    A 'gapGas' input file as for the case of gapFrapcon is always needed with the snapshot
    of gap compositions and pressure at time t = 0 or time of restart.

Usage
    In solverDict file:

    \verbatim
        gapGas timeTabulated;

        gapGasOptions
        {
            #include "$FOAM_CASE/constant/gapGasHistories";

            //timePoints (0.05000		604.42270	920.17348);
            //gapPressures (4.37114e+06	5.68517e+06	9.65693e+06);
            //Ar (0.0	0.0	0.0);
            //He (1.0	0.99959	0.87012);
            //Kr (0.0	0.00006	0.01987);
            //Ne (0.0	0.0	0.0);
            //Rn (0.0	0.0	0.0);
            //Xe (0.0	0.00035	0.11001);

            timeInterpolationMethod linear;
        }

    \endverbatim

    In gapGas file inside the time step folder (e.g. 0/ or 500.25/uniform/, etc):
    \verbatim
        massFractions
        {
            Ar 0;
            He 1;
            Kr 0;
            Ne 0;
            Rn 0;
            Xe 0;
        };

        gasPressure                     2.25e6;
    \endverbatim

    If using an external file like "gapGasHistories" dictionary, then in this file:

    \verbatim
        timePoints (0.0 100 1000 );
        gapPressures (4.37114e+06 5.68517e+06 9.65693e+06);
        Ar (0.0 0.0 0.0);
        He (1.0 0.99959 0.87012);
        Kr (0.0 0.00006 0.01987);
        Ne (0.0 0.0 0.0);
        Rn (0.0 0.0 0.0);
        Xe (0.0 0.00035 0.11001);
    \endverbatim

SourceFiles
    gapGasTimeTabulated.C

\mainauthor
Lokesh Verma, PSI

\contribution

\date 
    September 2024

\*---------------------------------------------------------------------------*/

#ifndef gapGasTimeTabulated_H
#define gapGasTimeTabulated_H

#include "gapGasModel.H"
#include "Function1.H"
#include "globalOptions.H"
#include "InterpolateTables.H"
#include <set>
#include "Table.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class gapFrapcon Declaration
\*---------------------------------------------------------------------------*/

class gapGasTimeTabulated
:
    public gapGasModel
{
    // Private typedefs

        typedef interpolateTableBase::interpolationMethod interpolationMethod;
        typedef interpolateTableBase::outOfBoundsMethod outOfBoundsMethod;

    // Private data

        //- Tabulated time-dependent gap gas pressure (Pa)
        scalarField pressureData_;

        //- Tabulated time-dependent mass fravtions for species
        scalarField ArgonData_;
        scalarField HeliumData_;
        scalarField KryptonData_;
        scalarField NeonData_;
        scalarField RadonData_;
        scalarField XenonData_;

        //- Tabulated time values
        scalarField timeData_;

        //- Method for time interpolation
        interpolationMethod timeMethod_;

        //- Method for out of bounds treatement
        outOfBoundsMethod boundsMethod_;

        //- Time-dependent table for gap gas pressure (Pa)
        scalarInterpolateTable pressureTable_;

        //- Time-dependent table for Argon mass fraction
        scalarInterpolateTable ArgonTable_;

        //- Time-dependent table for Helium mass fraction
        scalarInterpolateTable HeliumTable_;

        //- Time-dependent table for Krypton mass fraction
        scalarInterpolateTable KryptonTable_;

        //- Time-dependent table for Neon mass fraction
        scalarInterpolateTable NeonTable_;

        //- Time-dependent table for Radon mass fraction
        scalarInterpolateTable RadonTable_;

        //- Time-dependent table for Xenon mass fraction
        scalarInterpolateTable XenonTable_;

    // Private Member Functions

        //- Override the same functions in gapGasModel.
        virtual void correctMassFractions();

        virtual void correctMolFractions();

        //- Lookup the Pressures from the supplied physical time in seconds
        inline scalar lookupPressure(scalar t) const
        {
            return pressureTable_(mesh_.time().timeToUserTime(t));
        }

        //- Lookup the Argon mass fraction from the supplied physical time in seconds
        inline scalar lookupAr(scalar t) const
        {
            return ArgonTable_(mesh_.time().timeToUserTime(t));
        }

        //- Lookup the Helium mass fraction from the supplied physical time in seconds
        inline scalar lookupHe(scalar t) const
        {
            return HeliumTable_(mesh_.time().timeToUserTime(t));
        }

        //- Lookup the Krypton mass fraction from the supplied physical time in seconds
        inline scalar lookupKr(scalar t) const
        {
            return KryptonTable_(mesh_.time().timeToUserTime(t));
        }

        //- Lookup the Neon mass fraction from the supplied physical time in seconds
        inline scalar lookupNe(scalar t) const
        {
            return NeonTable_(mesh_.time().timeToUserTime(t));
        }

        //- Lookup the Radon mass fraction from the supplied physical time in seconds
        inline scalar lookupRn(scalar t) const
        {
            return RadonTable_(mesh_.time().timeToUserTime(t));
        }

        //- Lookup the Xenon mass fraction from the supplied physical time in seconds
        inline scalar lookupXe(scalar t) const
        {
            return XenonTable_(mesh_.time().timeToUserTime(t));
        }

        //- Disallow default bitwise copy construct
        gapGasTimeTabulated(const gapGasTimeTabulated&);

        //- Disallow default bitwise assignment
        void operator=(const gapGasTimeTabulated&);
        
protected:  

    // Protected data 

        //- Read As parameters for Frapcon conductivity model
        scalarList conductivityAs_;

        //- Read Bs parameters for Frapcon conductivity model
        scalarList conductivityBs_;

        //- Time dependent list of gas pressure values
        //autoPtr<Function1s::Table<scalar>> gasPressureList;
/*
        //- Time dependent list of gas composition values for gap gases
        autoPtr<Function1s::Table<scalar>> ArgonList;
        autoPtr<Function1s::Table<scalar>> HeliumList;
        autoPtr<Function1s::Table<scalar>> KryptonList;
        autoPtr<Function1s::Table<scalar>> NeonList;
        autoPtr<Function1s::Table<scalar>> RadonList;
        autoPtr<Function1s::Table<scalar>> XenonList;
*/
        //- gap gas species for mass fractions
        scalar Ar_;
        scalar He_;
        scalar Kr_;
        scalar Ne_;
        scalar Rn_;
        scalar Xe_;

public:           

    //- Runtime type information
    TypeName("timeTabulated");
    

    // Constructors

    //- Construct from mesh, materials and dict
    gapGasTimeTabulated
    (
        const fvMesh& mesh,
        const materials& mat,
        const volScalarField& T,
        const volVectorField& DorDD,
        const fissionGasRelease& fgr,
        const dictionary& gapGasDict
    );
        

    //- Destructor
    ~gapGasTimeTabulated();


    // Member Functions

    //- Gap gas thermal conductivity (from Frapcon manual)
        virtual scalar kappa(const scalar T) const;

    //- Return the accommodation coefficient (from Frapcon manual)
        virtual scalar a(const scalar T) const;

    //- Correct the gas thermodynamic properties
        virtual void correct();

            
    // IO
        
        //- Read data from time directory
        virtual bool readData(Istream& is);

        //- Read gapGasHistories from constant directory
        virtual void readHistories();

        //- Write data to time directory
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
