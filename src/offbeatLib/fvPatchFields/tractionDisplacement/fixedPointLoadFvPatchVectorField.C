/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "fixedPointLoadFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedPointLoadFvPatchVectorField::
fixedPointLoadFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(p, iF),
    point_(vector::zero),
    force_(vector::zero)
{}


fixedPointLoadFvPatchVectorField::
fixedPointLoadFvPatchVectorField
(
    const fixedPointLoadFvPatchVectorField& fplpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    tractionDisplacementFvPatchVectorField(fplpvf, p, iF, mapper),
    point_(fplpvf.point_),
    force_(fplpvf.force_)
{}


fixedPointLoadFvPatchVectorField::
fixedPointLoadFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    tractionDisplacementFvPatchVectorField(p, iF, dict, false),
    point_(dict.lookup("point")),
    force_(dict.lookup("force"))
{    
    if( mag(force_) == 0 )
    {
        FatalErrorInFunction
        << "The point load in " << patch().name() << " has zero magnitude." 
        << abort(FatalError);
    }

    // Read initial pressure value if present (for accurate restart)
    if(dict.found("pressure"))
    {
        pressure_ = scalarField("pressure", dict, p.size());
    }
}

fixedPointLoadFvPatchVectorField::
fixedPointLoadFvPatchVectorField
(
    const fixedPointLoadFvPatchVectorField& fplpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    tractionDisplacementFvPatchVectorField(fplpvf, iF), 
    point_(fplpvf.point_),
    force_(fplpvf.force_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedPointLoadFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    tractionDisplacementFvPatchVectorField::autoMap(m);
}


void fixedPointLoadFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    tractionDisplacementFvPatchVectorField::rmap(ptf, addr);
}


void fixedPointLoadFvPatchVectorField::updateTraction()
{
    vectorField n(patch().nf());

	const polyPatch& thisPatch = patch().patch();

	// List of patch points
	labelList labelPointsPatch( thisPatch.meshPoints() );

	// List of meshPoints
	pointField meshPoints(patch().boundaryMesh().mesh().points());

	scalar totalArea(0);

    // Find faces attached or containing the point given in point_. 
    // Distribute the input force over these faces as a stress.
	forAll(thisPatch, faceI)
	{
		face patchFace = thisPatch[faceI];

		forAll(patchFace , pointFaceI)
		{
			if(mag(meshPoints[patchFace[pointFaceI]] - point_) < 1e-6)
			{
				totalArea += patch().magSf()[faceI];
                vector forceDir = force_/mag(force_);
                scalar cosine = forceDir&n[faceI];
				pressure_[faceI] = -(mag(force_)/cosine);

				break;
			}
			else
            {
                pressure_[faceI] = 0;
            }
		}
	}

    reduce(totalArea, sumOp<scalar>());   
    
    pressure_ /= totalArea;
    traction_ = vector::zero;    
}


void fixedPointLoadFvPatchVectorField::write(Ostream& os) const
{
    tractionDisplacementFvPatchVectorField::write(os);

    os.writeKeyword("point") << point_
        << token::END_STATEMENT << nl;

    os.writeKeyword("force") << force_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fixedPointLoadFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
