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
    Foam::constantHeatSource

Description
	The `constant` heatSource class in OFFBEAT allows the user to set the volumetric
    heat source field in a file located the starting time folder. If the `Q` file is
	present in the starting time folder, the internal field and boundary conditions
	are set according to the `Q` file. Otherwise, the heat source field is set
	to 0 <code>W/m<sup>3</sup></code> in each cell, with zeroGradient BCs in all non-
	empty/-wedge patches. Note that when selecting this heat source model, the heat
	source field will remain the same as it was defined in the starting time folder
	(or as it was set by default by the code).

	\note
	At the moment, the constant model is the correct heat source model to use
	if you are coupling OFFBEAT to an external neutronics solver that provides the
	`Q` field directly. At every coupled iteration, simply substitute the
	`Q` file generated with the external solver in the starting time folder and
	(re)start the OFFBEAT simulation.

    \warning
    For backward-compatibility, `fromLatestTime` can be still used as the typename
    for the `constant` heatSource model. However, its use is deprecated.

________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for the `constant` heat
    source model

	\verbatim
        heatSource constant;
	\endverbatim

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

#ifndef constantHeatSource_H
#define constantHeatSource_H

#include "primitiveFields.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "materials.H"
#include "heatSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class constantHeatSource Declaration
\*---------------------------------------------------------------------------*/

class constantHeatSource
:
    public heatSource
{
    // Private data


    // Private Member Functions

        //- Disallow default bitwise copy construct
        constantHeatSource(const constantHeatSource&);

        //- Disallow default bitwise assignment
        void operator=(const constantHeatSource&);

 protected:

    // Protected data

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
    TypeName("constant");

    // Constructors

        //- Construct from mesh, materials and dict
        constantHeatSource
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceDict
        );


    // Selectors


    //- Destructor
    virtual ~constantHeatSource();


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
