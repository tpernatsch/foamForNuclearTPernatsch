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
    Foam::hydrideReorientationDesquines

Description:
    Hydride reorientation model using Desquines et al. [1]
    
    [1] J. Desquines, D. Drouan, M. Billone, M. Puls, P. March, S. Fourgeaud, C. Getrey, V. Elbaz, and M. Philippe,
    “Influence of temperature and hydrogen content on stress-induced radial hydride precipitation in Zircaloy-4 cladding”
    Journal of Nuclear Materials, vol. 453, no. 1, pp. 131–150, 2014

Usage
    In solverDict:
    \verbatim
    elementTransportOptions
    {
        solvers ( HNGD );

        HSolverOptions
        {
            stressDrivenDiffusion   true;
            hydrideReorientationModel   hydrideReorientationDesquines;

            materials
            {
                ...
            }
        }
    }
    \endverbatim

SourceFiles
    hydrideReorientationDesquines.C

L. Verma - Paul Scherrer Institut (PSI)

\contribution
    

\date 
    July 2024

\*---------------------------------------------------------------------------*/

#ifndef hydrideReorientationDesquines_H
#define hydrideReorientationDesquines_H

#include "hydrideReorientationModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class hydrideReorientationDesquines Declaration
\*---------------------------------------------------------------------------*/

class hydrideReorientationDesquines
:
    public hydrideReorientationModel
{
    // Private data
               
        //- Threshold stress value for radial reorientation
        volScalarField sigma0_;

        //- Stress level for 100% radial hydrides
        volScalarField sigma100_;

        //Fractional change in hydrides
        volScalarField fractionalHydridesChange_;

        //- Temperature required to dissolve Ctot during heatup
        volScalarField T_tssd_;
   
    // Private Member Functions
        
        //- Disallow default bitwise copy construct
        hydrideReorientationDesquines(const hydrideReorientationDesquines&);

        //- Disallow default bitwise assignment
        void operator=(const hydrideReorientationDesquines&);
     
protected:

    // Protected data

public:

    //- Runtime type information
    TypeName("hydrideReorientationDesquines");

    // Constructors

        //- Construct from mesh, materials and dictionary
        hydrideReorientationDesquines
        (
            const fvMesh& mesh,
            const dictionary& lawDict
        );

    // Selectors
    
    //- Destructor
    virtual ~hydrideReorientationDesquines();

    // Access Functions

    // Member Functions

    //- Correct
    virtual void correctReorientation
    (
        const labelList& addr
    );

    //- Calculate hydride reorientation
    //virtual void calculateReorientation(const scalar&, const scalar&, const scalar&);
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
