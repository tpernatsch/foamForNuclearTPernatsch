/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Description
    Central-differencing interpolation scheme class

\*---------------------------------------------------------------------------*/

#include "fvMesh.H"
#include "linearRegionCoupled.H"
#include "regionCoupledOFFBEATFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceInterpolationScheme(linearRegionCoupled)
}


template<class Type>
Foam::tmp<Foam::GeometricField<Foam::scalar, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::linearRegionCoupled<Type>::weights
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
  tmp<surfaceScalarField> tmpWeight
  (
    new surfaceScalarField
    (
      this->mesh().surfaceInterpolation::weights()
    )
  );

  forAll(tmpWeight->boundaryField(), patchI)
  {
      const fvPatch& patch = vf.boundaryField()[patchI].patch();
      scalarField& wP = tmpWeight->boundaryFieldRef()[patchI];
      
      if(isType<regionCoupledOFFBEATFvPatch>(patch))
      {  
        forAll(this->mesh().thisDb().names("volVectorField"), fieldI)
        {
          const word fieldName = this->mesh().thisDb().names("volVectorField")[fieldI];
          const fvPatchVectorField& fp = patch.lookupPatchField<volVectorField, vector>(fieldName);
          
          if(
              isType<regionCoupledOFFBEATFvPatchVectorField2>
              (
                fp
              )
              and
              !(
                fieldName.size() > 8 &&
                fieldName(fieldName.size()-8, 8) == "PrevIter"
              )
            )
          {
            const regionCoupledOFFBEATFvPatchVectorField2& rcPatchField = 
            refCast<const regionCoupledOFFBEATFvPatchVectorField2>(fp);

            wP = rcPatchField.giveWeights();

            Info << "wP is " << wP << endl;
          } 
        }
      }
  }    


            Info << "tmpWeight " << tmpWeight() << endl;
  return tmpWeight;
}

// ************************************************************************* //
