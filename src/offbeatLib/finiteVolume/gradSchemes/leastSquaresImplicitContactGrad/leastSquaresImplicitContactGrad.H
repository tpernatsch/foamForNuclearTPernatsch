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
    Foam::fv::leastSquaresImplicitContactGrad
 
Description
    Second-order gradient scheme using least-squares, specific for implicit
    contact BCs.
    
    Derived from solids4Foam and foamExtend.
 
SourceFiles
    leastSquaresImplicitContactGrad.C
 
 \*---------------------------------------------------------------------------*/
 
 #ifndef leastSquaresImplicitContactGrad_H
 #define leastSquaresImplicitContactGrad_H
 
 #include "gradScheme.H"
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace Foam
 {
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 namespace fv
 {
 
 /*---------------------------------------------------------------------------*\
                       Class leastSquaresImplicitContactGrad Declaration
 \*---------------------------------------------------------------------------*/
 
 template<class Type>
 class leastSquaresImplicitContactGrad
 :
     public fv::gradScheme<Type>
 {
     // Private Member Functions
 
         //- Disallow default bitwise copy construct
         leastSquaresImplicitContactGrad(const leastSquaresImplicitContactGrad&);
 
         //- Disallow default bitwise assignment
         void operator=(const leastSquaresImplicitContactGrad&);
 
 
 public:
 
     //- Runtime type information
     TypeName("leastSquaresImplicitContact");
 
 
     // Constructors
 
         //- Construct from mesh
         leastSquaresImplicitContactGrad(const fvMesh& mesh)
         :
             gradScheme<Type>(mesh)
         {}
 
         //- Construct from Istream
         leastSquaresImplicitContactGrad(const fvMesh& mesh, Istream&)
         :
             gradScheme<Type>(mesh)
         {}
 
 
     // Member Functions
 
         //- Return the gradient of the given field to the gradScheme::grad
         //  for optional caching
         virtual tmp
         <
             GeometricField
             <typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
         > calcGrad
         (
             const GeometricField<Type, fvPatchField, volMesh>& vsf,
             const word& name
         ) const;
 };
 
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace fv
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 } // End namespace Foam
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #ifdef NoRepository
     #include "leastSquaresImplicitContactGrad.C"
 #endif
 
 // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
 
 #endif
 
 // ************************************************************************* //