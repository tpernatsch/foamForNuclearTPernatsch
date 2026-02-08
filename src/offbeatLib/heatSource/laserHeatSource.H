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
    Foam::laserHeatSource

Description
    This heat source model is used to simulate the heating due to energy 
    deposition from an impinging laser on one of the surfaces of the domain. The
    model initialize a `volScalarField` for the laser intensity, called \f$I\f$.

    To compute the energy deposition in the domain, the model solves the 
    following radiation attenuation equation:

    $$
	\begin{aligned}
	\nabla \cdot (I \vec{n}) + \alpha I = 0
	\end{aligned}
	$$
	
	where:
	
	- \f$I\f$ is the radiation intensity \f$\left(m^{-2}s^{-1}\right)\f$;
	- \f$\vec{n}\f$ is the beam's direction unit vector;
	- \f$\alpha\f$ is the material absorption coefficient  \f$\left(m^{-1}\right)\f$.
	
    Once the model computes the spatial distribution of \f$I\f$, the resulting 
    power density distribution \f$Q\f$ as:

    $$
	\begin{aligned}
	Q = I \cdot \alpha
	\end{aligned}
	$$

    \f$\alpha\f$ is the absorption coefficient of the fully-dense material. 
    However, the code allows to take into account variation in local porosity
    of the material that can affect its absorption capabilities. This is done
    through an effective absorption coefficient \f$\alpha_{eff}\f$ that is 
    computed as:

    $$
	\begin{aligned}
	\alpha_{eff} = \alpha \cdot \frac{1-p}{1+0.5p}
	\end{aligned}
	$$

    where:

	- \f$p\f$ is the fractional local porosity of the material.

________________________________________________________________________________

\par Options
	The `laserHeatSource` heat source model requires the user to specify a few
	additional parameters in the `heatSourceOptions` sub-dictionary located 
	within the `solverDict`.
    
    Parameters in `heatSourceOptions`:
	
	- <b>`alphaAbsorption`</b> - The absorption coefficient of the material in \f$m^{-1}\f$ .
	- <b>`beamDirection`</b> - A vector that specifies the direction of the laser beam.
	- <b>`useEffectiveAlpha`</b> - A bool that specifies whether or not the effective 
    absorption coefficient should be used.
	- <b>`materials`</b> - A list of material (or cellZones names) where the heat source
	model applies.
	
________________________________________________________________________________

Usage
	Here is a code snippet of the `solverDict` to be used for the `laserHeatSource`
    heat source model:
	
	\verbatim
        heatSource laserHeatSource;

        heatSourceOptions
        {
            useEffectiveAlpha   on;
            alphaAbsorption     1e4;
            beamDirection       (0 0 -1);
            
            materials   (fuel);
        }
	\endverbatim
    
________________________________________________________________________________

SourceFiles
    laserHeatSource.C

\mainauthor
    E. Brunetto - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland
    Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    M. Reymond, A. Scolaro - EPFL

\date 
    November 2021

\*---------------------------------------------------------------------------*/

#ifndef laserHeatSource_H
#define laserHeatSource_H

#include "constantHeatSource.H"
#include "axialProfile.H"
#include "radialProfile.H"
#include "azimuthalProfile.H"
#include "volFields.H"
#include "InterpolateTables.H"
#include "globalOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class laserHeatSource Declaration
\*---------------------------------------------------------------------------*/

class laserHeatSource
:
    public constantHeatSource
{
 
    // Private data
    
    // Private Member Functions

        //- Disallow default bitwise copy construct
        laserHeatSource(const laserHeatSource&);

        //- Disallow default bitwise assignment
        void operator=(const laserHeatSource&);
        
protected:   

    // Protected data
        
        //- Laser Intensity (1/m2 s)
        volScalarField I_;

        //- Effective photon abs. cross section (1/m)
        scalar alphaAbs_;

        //- Effective photon abs. cross section Field(1/m)
        volScalarField alphaAbsField_;

        //- Switch to modify alpha with porosity
        Switch useEffectiveAlpha_;

        //- Beam direction (-)
        vector beamDirection_;

        //- Beam direction field (-)
        volVectorField beamDirectionField_;

        //- Addressing for power-producing cells
        labelList addr_;

public:

    //- Runtime type information
    TypeName("laserHeatSource");

    // Constructors

        //- Construct from mesh, materials and dict
        laserHeatSource
        (
            const fvMesh& mesh,
            const materials& materials,
            const dictionary& heatSourceDict,
            const bool timeDependent=false
        );


    //- Destructor
    ~laserHeatSource();


    // Member Functions

        //- Get addressing for power-producing cells
        void calcAddressing(const dictionary& dict);
        
        //- Calculate the fuel height
        void calcReferenceDimensions(const dictionary& dict);
    
        //- update the power distribution in the supplied list of cells
        virtual void correct();
        
        //- Return the next average power density (assuming the deltaT stays the
        // same)
        virtual scalar predictNextPowerDensity() const;    

        //- Return the average power density 
        virtual scalar averagePowerDensity() const;

        //- Correct effective absorption XS according to porosity
        virtual void correctEffectiveAbsorption();

        //- Reshape Q profile according to porosity
        virtual void reshapePowerPorosity();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
