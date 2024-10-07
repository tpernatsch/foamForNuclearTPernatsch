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
    Foam::plenumSpringFixedDisplacementFvPatchVectorField
 
Description
    Plenum spring patchfield for the top of the fuel column or for the topCap
    inner patch (fixedDisplacement version).

    The component of the displacement normal to the patch is derived from the
    force balance between internal stresses, spring force and plenum pressure.
    A plane strain (epsilonZ costant on the patch) or plane stress approximations
    (sigmaZ costant on the patch) are applied.
    The shear components of the displacement are set to zero.
    A pre-loading can be provided to the spring.

Usage
    \verbatim
    cladTop
    {
        type            plenumSpringFixedDisplacement;
        planeStrain     on;

        springModulus   1e4;
        initialSpringLoading 0;

        fuelTopPatches (fuelTop);
        topCapInnerPatches (cladOuter);

        value           $internalField;
    }    
    \endverbatim    
 
SourceFiles
    plenumSpringFixedDisplacementFvPatchVectorField.C

\mainauthor
    A. Scolaro - EPFL (ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, 
    Switzerland, Laboratory for Reactor Physics and Systems Behaviour)

\contribution
    E. Brunetto, C. Fiorina - EPFL\n
    I. Clifford - PSI (Paul Scherrer Institut, Switzerland)

\date 
    November 2021
 
\*---------------------------------------------------------------------------*/
 
#ifndef plenumSpringFixedDisplacementFvPatchVectorField_H
#define plenumSpringFixedDisplacementFvPatchVectorField_H
 
#include "fixedDisplacementFvPatchVectorField.H"
#include "plenumSpringBase.H"
 
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
namespace Foam
{
 
/*---------------------------------------------------------------------------*\
        Class plenumSpringFixedDisplacementFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/
 
class plenumSpringFixedDisplacementFvPatchVectorField
:
    public fixedDisplacementFvPatchVectorField,
    public plenumSpringBase
{
    // Private Data

protected:

    // Protected Data

        //- Apply plane strain approximation to deduce axial displacement
        Switch planeStrain_;


    // Protected Member Functions
 
public:
 
    //- Runtime type information
    TypeName("plenumSpringFixedDisplacement");


    // Constructors

        //- Construct from patch and internal field
        plenumSpringFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        plenumSpringFixedDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&,
            const bool valueRequired=true
        );

        //- Construct by mapping the given plenumSpringFixedDisplacementFvPatchVectorField
        //  onto a new patch
        plenumSpringFixedDisplacementFvPatchVectorField
        (
            const plenumSpringFixedDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&,
            const bool mappingRequired=true
        );


    // Member Functions

        //- Write
        virtual void write(Ostream&) const;
  
    // Evaluation functions

        //- Return the matrix diagonal coefficients corresponding to the
        //  evaluation of the value of this patchField with given weights
        virtual void updateCoeffs ();

 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //